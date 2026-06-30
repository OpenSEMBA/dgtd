#include "mor_rcs_driver.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
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

std::vector<std::pair<int, std::filesystem::path>> collectXFiles(const std::string& x_dir)
{
    std::vector<std::pair<int, std::filesystem::path>> files;
    const std::regex re(R"(^x_(\d+)$)");

    for (const auto& entry : std::filesystem::directory_iterator(x_dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        std::smatch m;
        const auto name = entry.path().filename().string();
        if (std::regex_match(name, m, re)) {
            files.emplace_back(std::stoi(m[1].str()), entry.path());
        }
    }

    std::sort(files.begin(), files.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    if (files.empty()) {
        throw std::runtime_error("No x_k files found in: " + x_dir);
    }

    return files;
}

std::pair<double, mfem::Vector> loadSnapshot(const std::filesystem::path& file_path,
                                             int expected_size)
{
    std::ifstream ifs(file_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open snapshot: " + file_path.string());
    }

    double time = 0.0;
    long long n = -1;
    if (!(ifs >> time)) {
        throw std::runtime_error("Invalid time header in: " + file_path.string());
    }
    if (!(ifs >> n)) {
        throw std::runtime_error("Invalid size header in: " + file_path.string());
    }

    if (n != expected_size) {
        throw std::runtime_error(
            "Snapshot size mismatch in " + file_path.string() +
            ". Expected " + std::to_string(expected_size) +
            ", got " + std::to_string(n));
    }

    mfem::Vector state(expected_size);
    for (int i = 0; i < expected_size; ++i) {
        if (!(ifs >> state[i])) {
            throw std::runtime_error(
                "Unexpected EOF while reading state values in: " + file_path.string());
        }
    }

    return {time, std::move(state)};
}

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

} // namespace

void runMORRCSPostProcessing(const MORRCSConfig& cfg)
{
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
    }

    const std::string data_path =
        cfg.output_root + "/RCSSurface/" + cfg.probe_name + "/";
    const std::string rank_path = data_path + "rank0";

    RCSSurfaceBinaryWriter writer(tags, &fec, fes, fields, rank_path);

    const int expected_size = fields.fieldBlockSize();
    const auto snapshots = collectXFiles(cfg.x_dir);

    int written = 0;
    for (const auto& [cycle, path] : snapshots) {
        (void)cycle;
        const auto [time, state] = loadSnapshot(path, expected_size);
        if (cfg.max_time.has_value() && time > cfg.max_time.value()) {
            continue;
        }
        fields.allDOFs() = state;
        writer.writeSnapshot(time);
        ++written;
    }

    if (written == 0) {
        throw std::runtime_error("No snapshots written (check --max-time or xdir contents).");
    }

    if (mfem::Mpi::WorldRank() == 0) {
        std::cout << "mor2rcs: wrote " << written << " surface snapshots to " << rank_path << '\n';
        std::cout << "mor2rcs: running RCS post-processing...\n";
    }

    std::vector<Frequency> frequencies = cfg.frequencies;
    RCSSurfacePostProcessor pp(data_path, cfg.case_path, frequencies, cfg.angles, cfg.max_time);
}

} // namespace maxwell::driver
