#include <mfem.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "components/FarField.h"
#include "driver/mor_rcs_driver.h"

namespace {

using json = nlohmann::json;

struct LauncherConfig {
    std::string case_path;
    std::string mesh_path;
    std::string x_dir;
    std::string ur_path;
    std::string xr_dir;
    std::string meta_path;
    std::string dump_x_dir;
    std::string output_root{"rcs_output"};
    std::string probe_name{"mor_rcs"};
    std::string rcs_json_path;
    std::vector<int> tags;
    std::optional<double> max_time;
    bool x_dir_explicit = false;
};

void printHelp()
{
    std::cout
        << "OpenSEMBA/dgtd mor2rcs\n"
        << "Replay MOR snapshots into rcssurface format and compute RCS.\n\n"
        << "Modes:\n"
        << "  Legacy:     --xdir <xfull/>          ASCII full-order x_k\n"
        << "  Reconstruct: --ur <Ur.bin> --xrdir <xr/>   x = Ur * xr on the fly\n\n"
        << "Usage:\n"
        << "  mor2rcs [--case <case.json>] [--mesh <mesh.msh>]\n"
        << "          [--xdir <xfull/> | --ur <Ur.bin> --xrdir <xr/>]\n"
        << "          [--meta <meta.json>] [--dump-xdir <dir>]\n"
        << "          [--tags <t1> <t2> ...] [--name <probe_name>] [--out <output_root>]\n"
        << "          [-i <rcs_sweep.json>] [--max-time <t>]\n\n"
        << "Defaults (when run in a post-processing folder):\n"
        << "  --case: unique .json file in current directory\n"
        << "  --mesh: unique .msh file in current directory\n"
        << "  --xdir: ./xfull if present, else ./x  (legacy mode only)\n"
        << "  --meta: <dir(Ur)>/meta.json if present\n"
        << "  --name: mor_rcs\n"
        << "  --out : ./rcs_output\n\n"
        << "Ur.bin: float64, column-major (N×r), host endian. Peak RAM ≈ Ur + one x.\n";
}

std::string findUniqueFileWithExtension(const std::filesystem::path& dir,
                                        const std::string& ext)
{
    std::vector<std::filesystem::path> matches;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        if (entry.path().extension() == ext) {
            matches.push_back(entry.path());
        }
    }

    if (matches.empty()) {
        throw std::runtime_error("No file with extension '" + ext +
                                 "' found in current directory.");
    }
    if (matches.size() > 1) {
        throw std::runtime_error(
            "Multiple files with extension '" + ext +
            "' found in current directory. Please pass explicit path via CLI.");
    }
    return matches.front().string();
}

std::string resolveXDir(const std::filesystem::path& cwd)
{
    if (std::filesystem::is_directory(cwd / "xfull")) {
        return (cwd / "xfull").string();
    }
    if (std::filesystem::is_directory(cwd / "x")) {
        return (cwd / "x").string();
    }
    return (cwd / "xfull").string();
}

LauncherConfig parseArgs(int argc, char** argv)
{
    LauncherConfig cfg;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printHelp();
            std::exit(0);
        }
        else if (arg == "--case" && i + 1 < argc) {
            cfg.case_path = argv[++i];
        }
        else if (arg == "--mesh" && i + 1 < argc) {
            cfg.mesh_path = argv[++i];
        }
        else if (arg == "--xdir" && i + 1 < argc) {
            cfg.x_dir = argv[++i];
            cfg.x_dir_explicit = true;
        }
        else if (arg == "--ur" && i + 1 < argc) {
            cfg.ur_path = argv[++i];
        }
        else if (arg == "--xrdir" && i + 1 < argc) {
            cfg.xr_dir = argv[++i];
        }
        else if (arg == "--meta" && i + 1 < argc) {
            cfg.meta_path = argv[++i];
        }
        else if (arg == "--dump-xdir" && i + 1 < argc) {
            cfg.dump_x_dir = argv[++i];
        }
        else if (arg == "--name" && i + 1 < argc) {
            cfg.probe_name = argv[++i];
        }
        else if (arg == "--out" && i + 1 < argc) {
            cfg.output_root = argv[++i];
        }
        else if (arg == "-i" && i + 1 < argc) {
            cfg.rcs_json_path = argv[++i];
        }
        else if (arg == "--max-time" && i + 1 < argc) {
            cfg.max_time = std::stod(argv[++i]);
        }
        else if (arg == "--tags") {
            while (i + 1 < argc) {
                const std::string next = argv[i + 1];
                if (!next.empty() && next[0] == '-') {
                    break;
                }
                cfg.tags.push_back(std::stoi(argv[++i]));
            }
        }
        else {
            throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
    }

    const auto cwd = std::filesystem::current_path();
    if (cfg.case_path.empty()) {
        cfg.case_path = findUniqueFileWithExtension(cwd, ".json");
    }
    if (cfg.mesh_path.empty()) {
        cfg.mesh_path = findUniqueFileWithExtension(cwd, ".msh");
    }

    const bool ur_mode = !cfg.ur_path.empty() || !cfg.xr_dir.empty();
    if (ur_mode) {
        if (cfg.ur_path.empty() || cfg.xr_dir.empty()) {
            throw std::runtime_error(
                "Reconstruct mode requires both --ur <Ur.bin> and --xrdir <xr/>.");
        }
        if (cfg.x_dir_explicit) {
            std::cerr << "mor2rcs: warning: --xdir ignored because --ur/--xrdir is set\n";
            cfg.x_dir.clear();
        }
    } else if (cfg.x_dir.empty()) {
        cfg.x_dir = resolveXDir(cwd);
    }

    if (!std::filesystem::is_regular_file(cfg.case_path)) {
        throw std::runtime_error("Case JSON not found: " + cfg.case_path);
    }
    if (!std::filesystem::is_regular_file(cfg.mesh_path)) {
        throw std::runtime_error("Mesh file not found: " + cfg.mesh_path);
    }
    if (!ur_mode && !std::filesystem::is_directory(cfg.x_dir)) {
        throw std::runtime_error("x directory not found: " + cfg.x_dir);
    }
    if (ur_mode) {
        if (!std::filesystem::is_regular_file(cfg.ur_path)) {
            throw std::runtime_error("Ur file not found: " + cfg.ur_path);
        }
        if (!std::filesystem::is_directory(cfg.xr_dir)) {
            throw std::runtime_error("xr directory not found: " + cfg.xr_dir);
        }
        if (!cfg.meta_path.empty() && !std::filesystem::is_regular_file(cfg.meta_path)) {
            throw std::runtime_error("meta JSON not found: " + cfg.meta_path);
        }
    }
    if (cfg.rcs_json_path.empty()) {
        throw std::runtime_error("RCS sweep JSON required (-i <rcs_sweep.json>).");
    }

    return cfg;
}

std::vector<double> linspace(double a, double b, size_t N)
{
    if (N == 0) return {};
    if (N == 1) return {a};
    std::vector<double> xs(N);
    const double h = (b - a) / static_cast<double>(N - 1);
    for (size_t i = 0; i < N; ++i) {
        xs[i] = a + i * h;
    }
    return xs;
}

maxwell::driver::MORRCSConfig buildDriverConfig(const LauncherConfig& launcher_cfg)
{
    std::ifstream f(launcher_cfg.rcs_json_path);
    if (!f) {
        throw std::runtime_error("Cannot open RCS sweep JSON: " + launcher_cfg.rcs_json_path);
    }
    const auto rcs_input = json::parse(f);

    maxwell::driver::MORRCSConfig cfg;
    cfg.case_path = launcher_cfg.case_path;
    cfg.mesh_path = launcher_cfg.mesh_path;
    cfg.x_dir = launcher_cfg.x_dir;
    cfg.ur_path = launcher_cfg.ur_path;
    cfg.xr_dir = launcher_cfg.xr_dir;
    cfg.meta_path = launcher_cfg.meta_path;
    cfg.dump_x_dir = launcher_cfg.dump_x_dir;
    cfg.output_root = launcher_cfg.output_root;
    cfg.probe_name = launcher_cfg.probe_name;
    cfg.tags = launcher_cfg.tags;
    cfg.max_time = launcher_cfg.max_time;

    const auto& freq_spec = rcs_input.at("frequencies");
    cfg.frequencies = linspace(
        freq_spec.at("start").get<double>(),
        freq_spec.at("end").get<double>(),
        freq_spec.at("steps").get<size_t>());

    const auto& ang_spec = rcs_input.at("angles");
    const auto thetas = linspace(
        ang_spec.at("theta").at("start").get<double>(),
        ang_spec.at("theta").at("end").get<double>(),
        ang_spec.at("theta").at("steps").get<size_t>());
    const auto phis = linspace(
        ang_spec.at("phi").at("start").get<double>(),
        ang_spec.at("phi").at("end").get<double>(),
        ang_spec.at("phi").at("steps").get<size_t>());

    for (const auto t : thetas) {
        for (const auto p : phis) {
            cfg.angles.push_back({t, p});
        }
    }

    if (cfg.max_time.has_value() == false &&
        rcs_input.contains("max_time") && !rcs_input["max_time"].is_null()) {
        cfg.max_time = rcs_input["max_time"].get<double>();
    }

    return cfg;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        mfem::Mpi::Init(argc, argv);
        mfem::Hypre::Init();

        if (mfem::Mpi::WorldSize() != 1) {
            throw std::runtime_error("mor2rcs only supports single-rank execution.");
        }

        mfem::Device device("cpu");
        if (mfem::Mpi::WorldRank() == 0) {
            device.Print();
        }

        const LauncherConfig launcher_cfg = parseArgs(argc, argv);
        const auto driver_cfg = buildDriverConfig(launcher_cfg);

        maxwell::driver::runMORRCSPostProcessing(driver_cfg);

        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "mor2rcs completed. Output in '" << driver_cfg.output_root << "'.\n";
        }

        mfem::Mpi::Finalize();
        return 0;
    }
    catch (const std::exception& e) {
        if (mfem::Mpi::WorldRank() == 0) {
            std::cerr << "mor2rcs error: " << e.what() << std::endl;
        }
        mfem::Mpi::Finalize();
        return 1;
    }
}
