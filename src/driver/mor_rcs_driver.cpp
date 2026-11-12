#include "mor_rcs_driver.h"

#include "MorUrReconstructor.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

#ifdef SEMBA_DGTD_ENABLE_OPENMP
#include <omp.h>
#endif

#include <nlohmann/json.hpp>

#include "components/RCSSurfaceBinaryWriter.h"
#include "components/RCSSurfacePostProcessor.h"
#include "driver/driver.h"
#include "evolution/Fields.h"

namespace maxwell::driver {

namespace {

using json = nlohmann::json;

std::vector<int> tagsFromRCSProbes(const json& case_data)
{
    std::unordered_set<int> unique_tags;
    if (!case_data.contains("probes") || !case_data["probes"].contains("rcssurface")) {
        return {};
    }

    for (const auto& probe : case_data["probes"]["rcssurface"]) {
        for (const auto& tag : probe.at("tags")) {
            unique_tags.insert(tag.get<int>());
        }
    }

    return std::vector<int>(unique_tags.begin(), unique_tags.end());
}

std::vector<int> tagsFromPlaneWaveSources(const json& case_data)
{
    if (!case_data.contains("sources")) {
        return {};
    }

    std::vector<std::vector<int>> source_tag_sets;
    for (const auto& src : case_data["sources"]) {
        if (!src.contains("type") || src["type"].get<std::string>() != "planewave") {
            continue;
        }
        if (!src.contains("tags")) {
            continue;
        }
        std::vector<int> tags;
        for (const auto& tag : src["tags"]) {
            tags.push_back(tag.get<int>());
        }
        std::sort(tags.begin(), tags.end());
        source_tag_sets.push_back(std::move(tags));
    }

    if (source_tag_sets.empty()) {
        return {};
    }

    const auto& first = source_tag_sets.front();
    for (size_t i = 1; i < source_tag_sets.size(); ++i) {
        if (source_tag_sets[i] != first) {
            throw std::runtime_error(
                "Multiple planewave sources with different tags; pass --tags explicitly.");
        }
    }

    return first;
}

std::vector<int> resolveSurfaceTags(const json& case_data, const std::vector<int>& cli_tags)
{
    if (!cli_tags.empty()) {
        return cli_tags;
    }

    auto tags = tagsFromRCSProbes(case_data);
    if (!tags.empty()) {
        return tags;
    }

    tags = tagsFromPlaneWaveSources(case_data);
    if (!tags.empty()) {
        return tags;
    }

    throw std::runtime_error(
        "Could not resolve NTFF surface tags. Add rcssurface probe tags, planewave source "
        "tags, or pass --tags on the command line.");
}

void dumpAsciiSnapshot(const std::string& dir,
                       int index,
                       double time,
                       const mfem::Vector& state)
{
    std::filesystem::create_directories(dir);
    const auto path = std::filesystem::path(dir) / ("x_" + std::to_string(index));
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("Cannot write debug snapshot: " + path.string());
    }
    ofs << std::scientific << std::setprecision(16);
    ofs << time << "\n" << state.Size() << "\n";
    for (int i = 0; i < state.Size(); ++i) {
        ofs << state[i] << "\n";
    }
}

int writeSnapshotsFromXDir(const MORRCSConfig& cfg,
                           maxwell::Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields,
                           RCSSurfaceBinaryWriter& writer,
                           int expected_size)
{
    const auto snapshots = collectIndexedSnapshots(cfg.x_dir, "x");
    int written = 0;
    for (const auto& [cycle, path] : snapshots) {
        const auto [time, state] = loadAsciiSnapshot(path, expected_size);
        if (cfg.max_time.has_value() && time > cfg.max_time.value()) {
            continue;
        }
        fields.allDOFs() = state;
        writer.writeSnapshot(time);
        if (!cfg.dump_x_dir.empty()) {
            dumpAsciiSnapshot(cfg.dump_x_dir, cycle, time, state);
        }
        ++written;
    }
    return written;
}

int writeSurfaceFromUrXr(
    std::vector<double>& Ur,
    const MorUrMeta& meta,
    const std::vector<std::pair<int, std::string>>& snapshots,
    const MORRCSConfig& cfg,
    maxwell::Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields,
    RCSSurfaceBinaryWriter& writer,
    int expected_size);

int writeSnapshotsFromUrXr(const MORRCSConfig& cfg,
                           maxwell::Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields,
                           RCSSurfaceBinaryWriter& writer,
                           int expected_size)
{
    std::string meta_path = cfg.meta_path;
    if (meta_path.empty()) {
        const auto sibling = std::filesystem::path(cfg.ur_path).parent_path() / "meta.json";
        if (std::filesystem::is_regular_file(sibling)) {
            meta_path = sibling.string();
        }
    }

    MorUrMeta meta;
    meta.N = expected_size;
    meta.layout = "colmajor";
    meta.dtype = "float64";

    if (!meta_path.empty()) {
        meta = loadMorUrMeta(meta_path);
        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "mor2rcs: loaded meta " << meta_path
                      << " (N=" << meta.N << ", r=" << meta.r
                      << ", layout=" << meta.layout << ")\n";
        }
    } else {
        const auto probe = collectIndexedSnapshots(cfg.xr_dir, "xr");
        std::ifstream ifs(probe.front().second);
        double t = 0.0;
        long long rpeek = -1;
        if (!(ifs >> t >> rpeek) || rpeek <= 0) {
            throw std::runtime_error(
                "Cannot infer r from " + probe.front().second +
                "; provide meta.json with N and r.");
        }
        meta.r = rpeek;
        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "mor2rcs: no meta.json; using N=" << meta.N
                      << " from FE space, r=" << meta.r << " from first xr\n";
        }
    }

    if (meta.N != expected_size) {
        throw std::runtime_error(
            "Ur rows N=" + std::to_string(meta.N) +
            " != fields.fieldBlockSize()=" + std::to_string(expected_size));
    }
    if (meta.r <= 0) {
        throw std::runtime_error("Invalid reduced rank r in Ur/meta.");
    }

    if (mfem::Mpi::WorldRank() == 0) {
        const double ur_gib =
            static_cast<double>(meta.N) * static_cast<double>(meta.r) * 8.0 /
            (1024.0 * 1024.0 * 1024.0);
        std::cout << "mor2rcs: loading Ur from " << cfg.ur_path
                  << " (~" << std::fixed << std::setprecision(2) << ur_gib
                  << " GiB for Ur alone; xfull not required)\n";
    }

    auto Ur = loadUrBinary(cfg.ur_path, meta.N, meta.r, meta.layout, meta.dtype);
    const auto snapshots = collectIndexedSnapshots(cfg.xr_dir, "xr");

    if (!cfg.dump_x_dir.empty()) {
        mfem::Vector x_full(expected_size);
        int written = 0;
        for (const auto& [cycle, path] : snapshots) {
            const auto [time, xr] = loadAsciiSnapshot(path, static_cast<int>(meta.r));
            if (cfg.max_time.has_value() && time > cfg.max_time.value()) {
                continue;
            }
            urMatVec(Ur, meta.N, meta.r, xr, x_full);
            fields.allDOFs() = x_full;
            writer.writeSnapshot(time);
            dumpAsciiSnapshot(cfg.dump_x_dir, cycle, time, x_full);
            ++written;
        }
        return written;
    }

    return writeSurfaceFromUrXr(Ur, meta, snapshots, cfg, fields, writer, expected_size);
}

void loadXrColumn(const std::string& path, int r, double* column, double& time)
{
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open snapshot: " + path);
    }
    long long n = -1;
    if (!(ifs >> time) || !(ifs >> n)) {
        throw std::runtime_error("Invalid xr header in: " + path);
    }
    if (n != r) {
        throw std::runtime_error(
            "Snapshot size mismatch in " + path + ". Expected " +
            std::to_string(r) + ", got " + std::to_string(n));
    }
    for (int j = 0; j < r; ++j) {
        if (!(ifs >> column[j])) {
            throw std::runtime_error("Unexpected EOF in: " + path);
        }
    }
}

int writeSurfaceFromUrXr(
    std::vector<double>& Ur,
    const MorUrMeta& meta,
    const std::vector<std::pair<int, std::string>>& snapshots,
    const MORRCSConfig& cfg,
    maxwell::Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields,
    RCSSurfaceBinaryWriter& writer,
    int expected_size)
{
    mfem::Array<int> ids;
    std::vector<double> signs;
    writer.surfaceParentMap(ids, signs);

    const int m = ids.Size();
    const int ndofs = fields.ndofsPerComponent();
    const int ncomp = 6;
    const int rows = ncomp * m;
    const int r = static_cast<int>(meta.r);
    const int nsnap = static_cast<int>(snapshots.size());

    if (m != writer.numSurfaceDofs()) {
        throw std::runtime_error(
            "Surface DOF map size " + std::to_string(m) +
            " != writer surface size " + std::to_string(writer.numSurfaceDofs()));
    }
    if (ncomp * ndofs != expected_size || meta.N != expected_size) {
        throw std::runtime_error("Field block size does not match Ur rows.");
    }
    for (int i = 0; i < m; ++i) {
        if (ids[i] < 0 || ids[i] >= ndofs) {
            throw std::runtime_error(
                "Surface parent DOF " + std::to_string(ids[i]) + " is outside 0.." +
                std::to_string(ndofs - 1));
        }
    }

    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "mor2rcs: RCS uses " << rows
                  << " surface unknowns, not the full " << expected_size << '\n';
    }

    std::vector<double> A(static_cast<std::size_t>(rows) * static_cast<std::size_t>(r));
#ifdef SEMBA_DGTD_ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < rows; ++row) {
        const int comp = row / m;
        const int i = row - comp * m;
        const std::int64_t g =
            static_cast<std::int64_t>(ids[i]) +
            static_cast<std::int64_t>(comp) * ndofs;
        const double s = signs[static_cast<std::size_t>(i)];
        double* dst = A.data() + static_cast<std::size_t>(row) * static_cast<std::size_t>(r);
        for (int j = 0; j < r; ++j) {
            dst[j] = s * Ur[static_cast<std::size_t>(g + static_cast<std::int64_t>(j) * meta.N)];
        }
    }

    {
        const int candidates[] = {
            std::max(0, nsnap / 5),
            std::max(0, nsnap / 2),
            std::max(0, (4 * nsnap) / 5)
        };
        double strongest_ref = 0.0;
        for (const int check_k : candidates) {
            const auto [time, xr] = loadAsciiSnapshot(
                snapshots[static_cast<std::size_t>(check_k)].second, r);
            mfem::Vector x_full;
            urMatVec(Ur, meta.N, r, xr, x_full);
            fields.allDOFs() = x_full;
            std::vector<double> transferred(static_cast<std::size_t>(rows));
            writer.copyTransferredSurface(transferred.data());
            double max_abs = 0.0;
            double max_ref = 0.0;
            for (int row = 0; row < rows; ++row) {
                const double* a = A.data() + static_cast<std::size_t>(row) * static_cast<std::size_t>(r);
                double sum = 0.0;
                for (int j = 0; j < r; ++j) {
                    sum += a[j] * xr[j];
                }
                const double ref = transferred[static_cast<std::size_t>(row)];
                max_abs = std::max(max_abs, std::abs(sum - ref));
                max_ref = std::max(max_ref, std::abs(ref));
            }
            strongest_ref = std::max(strongest_ref, max_ref);
            if (mfem::Mpi::WorldRank() == 0) {
                std::cout << std::scientific << std::setprecision(3)
                          << "mor2rcs: transfer check at snapshot " << check_k
                          << " max_abs=" << max_abs << " max_ref=" << max_ref << '\n';
            }
            if (max_abs > 1e-6 * std::max(1.0, max_ref)) {
                throw std::runtime_error(
                    "Surface Ur*xr does not match the RCS transfer map.");
            }
        }
        if (strongest_ref < 1e-8) {
            throw std::runtime_error(
                "Surface transfer check never saw a nonzero field.");
        }
    }

    std::vector<double>().swap(Ur);

    std::vector<double> times(static_cast<std::size_t>(nsnap));
    std::vector<double> X(static_cast<std::size_t>(r) * static_cast<std::size_t>(nsnap));
    std::atomic<int> bad{-1};
    std::string bad_msg;
    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "mor2rcs: reading " << nsnap << " reduced states\n";
    }
#ifdef SEMBA_DGTD_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic, 32)
#endif
    for (int k = 0; k < nsnap; ++k) {
        try {
            loadXrColumn(snapshots[static_cast<std::size_t>(k)].second,
                         r,
                         X.data() + static_cast<std::size_t>(k) * static_cast<std::size_t>(r),
                         times[static_cast<std::size_t>(k)]);
        } catch (const std::exception& ex) {
#ifdef SEMBA_DGTD_ENABLE_OPENMP
#pragma omp critical
#endif
            if (bad.load() < 0) {
                bad.store(k);
                bad_msg = ex.what();
            }
        }
    }
    if (bad.load() >= 0) {
        throw std::runtime_error(bad_msg);
    }

    std::vector<int> used;
    used.reserve(static_cast<std::size_t>(nsnap));
    for (int k = 0; k < nsnap; ++k) {
        if (cfg.max_time.has_value() && times[static_cast<std::size_t>(k)] > cfg.max_time.value()) {
            continue;
        }
        used.push_back(k);
    }

    const int batch = 512;
    std::vector<double> Y(static_cast<std::size_t>(rows) * static_cast<std::size_t>(batch));
    const auto t0 = std::chrono::steady_clock::now();
    int written = 0;
    for (int b0 = 0; b0 < static_cast<int>(used.size()); b0 += batch) {
        const int nb = std::min(batch, static_cast<int>(used.size()) - b0);
#ifdef SEMBA_DGTD_ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int row = 0; row < rows; ++row) {
            const double* a = A.data() + static_cast<std::size_t>(row) * static_cast<std::size_t>(r);
            for (int t = 0; t < nb; ++t) {
                const int k = used[static_cast<std::size_t>(b0 + t)];
                const double* x = X.data() + static_cast<std::size_t>(k) * static_cast<std::size_t>(r);
                double sum = 0.0;
                for (int j = 0; j < r; ++j) {
                    sum += a[j] * x[j];
                }
                Y[static_cast<std::size_t>(row) + static_cast<std::size_t>(t) * static_cast<std::size_t>(rows)] = sum;
            }
        }
        for (int t = 0; t < nb; ++t) {
            const int k = used[static_cast<std::size_t>(b0 + t)];
            writer.writePackedSnapshot(
                times[static_cast<std::size_t>(k)],
                Y.data() + static_cast<std::size_t>(t) * static_cast<std::size_t>(rows));
        }
        written += nb;
        if (mfem::Mpi::WorldRank() == 0 && (written == nb || written % 8192 < batch || written == static_cast<int>(used.size()))) {
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "mor2rcs: wrote " << written << " / " << used.size()
                      << " surface snapshots in " << std::fixed << std::setprecision(1)
                      << sec << " s\n";
        }
    }
    writer.flush();
    return written;
}

} // namespace

void runMORRCSPostProcessing(const MORRCSConfig& cfg)
{
    if (cfg.useUrReconstruct()) {
        if (!std::filesystem::is_regular_file(cfg.ur_path)) {
            throw std::runtime_error("Ur file not found: " + cfg.ur_path);
        }
        if (!std::filesystem::is_directory(cfg.xr_dir)) {
            throw std::runtime_error("xr directory not found: " + cfg.xr_dir);
        }
    } else {
        if (cfg.x_dir.empty() || !std::filesystem::is_directory(cfg.x_dir)) {
            throw std::runtime_error(
                "Provide --xdir <xfull/> or both --ur <Ur.bin> and --xrdir <xr/>.");
        }
    }

    auto case_data = parseJSONfile(cfg.case_path);
    case_data["model"]["filename"] = std::filesystem::path(cfg.mesh_path).filename().string();

    const auto mesh_path = std::filesystem::path(cfg.mesh_path);
    const std::string synthetic_case_path =
        (mesh_path.parent_path() / (mesh_path.stem().string() + ".json")).string();

    maxwell::Model model = buildModel(case_data, synthetic_case_path, false);
    maxwell::SolverOptions opts = buildSolverOptions(case_data);

    mfem::DG_FECollection fec(opts.evolution.order,
                              model.getMesh().Dimension(),
                              opts.basis_type);
    mfem::ParFiniteElementSpace fes(&model.getMesh(), &fec);
    fes.ExchangeFaceNbrData();
    fes.GetParMesh()->ExchangeFaceNbrData();

    maxwell::Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction> fields(fes);

    const auto tags = resolveSurfaceTags(case_data, cfg.tags);
    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "mor2rcs: NTFF surface tags:";
        for (const int tag : tags) {
            std::cout << ' ' << tag;
        }
        std::cout << '\n';
        if (cfg.useUrReconstruct() && cfg.dump_x_dir.empty()) {
            std::cout << "mor2rcs: mode=surface rows of Ur*xr\n";
        } else if (cfg.useUrReconstruct()) {
            std::cout << "mor2rcs: mode=Ur*xr reconstruct-on-the-fly\n";
        } else {
            std::cout << "mor2rcs: mode=legacy ASCII --xdir\n";
        }
    }

    const std::string data_path =
        cfg.output_root + "/RCSSurface/" + cfg.probe_name + "/";
    const std::string rank_path = data_path + "rank0";

    RCSSurfaceBinaryWriter writer(tags, &fec, fes, fields, rank_path);

    const int expected_size = fields.fieldBlockSize();
    const int written = cfg.useUrReconstruct()
        ? writeSnapshotsFromUrXr(cfg, fields, writer, expected_size)
        : writeSnapshotsFromXDir(cfg, fields, writer, expected_size);

    if (written == 0) {
        throw std::runtime_error(
            "No snapshots written (check --max-time or snapshot directory contents).");
    }

    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "mor2rcs: wrote " << written << " surface snapshots to " << rank_path << '\n';
        std::cout << "mor2rcs: running RCS post-processing...\n";
    }

    std::vector<Frequency> frequencies = cfg.frequencies;
    RCSSurfacePostProcessor pp(data_path, cfg.case_path, frequencies, cfg.angles, cfg.max_time);
}

} // namespace maxwell::driver
