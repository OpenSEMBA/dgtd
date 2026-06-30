#include <mfem.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "driver/driver.h"
#include "evolution/Fields.h"

namespace {

struct ReplayConfig {
    std::string case_path;
    std::string mesh_path;
    std::string x_dir;
    std::string output_root{"output"};
    std::string dataset_name{"mor_paraview.vtk"};
};

void printHelp()
{
    std::cout
        << "OpenSEMBA/dgtd mor2paraview\n"
        << "Replay MOR state vectors x_k into ParaView time-series output.\n\n"
        << "Usage:\n"
        << "  mor2paraview [--case <case.json>] [--mesh <mesh.msh>] [--xdir <x_folder>]\\n"
        << "          [--name <dataset_name>] [--out <output_folder>]\n\n"
        << "Defaults (when run in a post-processing folder):\n"
        << "  --case: unique .json file in current directory\n"
        << "  --mesh: unique .msh file in current directory\n"
        << "  --xdir: ./x\n"
        << "  --name: mor_paraview.vtk\n"
        << "  --out : ./output\n";
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

ReplayConfig parseArgs(int argc, char** argv)
{
    ReplayConfig cfg;

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
        }
        else if (arg == "--name" && i + 1 < argc) {
            cfg.dataset_name = argv[++i];
        }
        else if (arg == "--out" && i + 1 < argc) {
            cfg.output_root = argv[++i];
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
    if (cfg.x_dir.empty()) {
        cfg.x_dir = (cwd / "x").string();
    }

    if (!std::filesystem::is_regular_file(cfg.case_path)) {
        throw std::runtime_error("Case JSON not found: " + cfg.case_path);
    }
    if (!std::filesystem::is_regular_file(cfg.mesh_path)) {
        throw std::runtime_error("Mesh file not found: " + cfg.mesh_path);
    }
    if (!std::filesystem::is_directory(cfg.x_dir)) {
        throw std::runtime_error("x directory not found: " + cfg.x_dir);
    }

    return cfg;
}

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

} // namespace

int main(int argc, char** argv)
{
    try {
        mfem::Mpi::Init(argc, argv);
        mfem::Hypre::Init();

        if (mfem::Mpi::WorldSize() != 1) {
            throw std::runtime_error("mor2paraview only supports single-rank execution.");
        }

        mfem::Device device("cpu");
        if (mfem::Mpi::WorldRank() == 0) {
            device.Print();
        }

        const ReplayConfig cfg = parseArgs(argc, argv);

        auto case_data = maxwell::driver::parseJSONfile(cfg.case_path);
        case_data["model"]["filename"] = std::filesystem::path(cfg.mesh_path).filename().string();

        const auto mesh_path = std::filesystem::path(cfg.mesh_path);
        const std::string synthetic_case_path =
            (mesh_path.parent_path() / (mesh_path.stem().string() + ".json")).string();

        maxwell::Model model = maxwell::driver::buildModel(case_data, synthetic_case_path, false);
        maxwell::SolverOptions opts = maxwell::driver::buildSolverOptions(case_data);

        mfem::DG_FECollection fec(opts.evolution.order,
                                  model.getMesh().Dimension(),
                                  opts.basis_type);
        mfem::ParFiniteElementSpace fes(&model.getMesh(), &fec);
        fes.ExchangeFaceNbrData();
        fes.GetParMesh()->ExchangeFaceNbrData();

        maxwell::Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction> fields(fes);

        std::filesystem::create_directories(cfg.output_root);

        mfem::ParaViewDataCollection pd(cfg.dataset_name, fes.GetParMesh());
        pd.SetPrefixPath(cfg.output_root);
        pd.RegisterField("E", &fields.get(maxwell::FieldType::E));
        pd.RegisterField("H", &fields.get(maxwell::FieldType::H));

        const auto geom_order = fes.GetMesh()->GetElementTransformation(0)->Order();
        const auto fec_order = fes.FEColl()->GetOrder();
        const bool high_order = (geom_order > 1) || (fec_order > 1);
        pd.SetHighOrderOutput(high_order);
        pd.SetLevelsOfDetail(std::max(geom_order, fec_order));
        pd.SetDataFormat(mfem::VTKFormat::BINARY);

        const int expected_size = fields.fieldBlockSize();
        const auto snapshots = collectXFiles(cfg.x_dir);

        int written = 0;
        for (const auto& [cycle, path] : snapshots) {
            const auto [time, state] = loadSnapshot(path, expected_size);
            fields.allDOFs() = state;
            pd.SetCycle(cycle);
            pd.SetTime(time);
            pd.Save();
            written++;
        }

        if (mfem::Mpi::WorldRank() == 0) {
            std::cout << "mor2paraview completed. Wrote " << written
                      << " snapshots to '" << cfg.output_root << "'.\n";
        }

        mfem::Mpi::Finalize();
        return 0;
    }
    catch (const std::exception& e) {
        if (mfem::Mpi::WorldRank() == 0) {
            std::cerr << "mor2paraview error: " << e.what() << std::endl;
        }
        mfem::Mpi::Finalize();
        return 1;
    }
}
