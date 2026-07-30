#include "mor_rcs_driver.h"

#include "MorUrReconstructor.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

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

    const auto Ur = loadUrBinary(cfg.ur_path, meta.N, meta.r, meta.layout, meta.dtype);
    const auto snapshots = collectIndexedSnapshots(cfg.xr_dir, "xr");

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
        if (!cfg.dump_x_dir.empty()) {
            dumpAsciiSnapshot(cfg.dump_x_dir, cycle, time, x_full);
        }
        ++written;
    }
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
        if (cfg.useUrReconstruct()) {
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
