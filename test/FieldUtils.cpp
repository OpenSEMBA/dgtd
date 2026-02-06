#include "FieldUtils.h"
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace maxwell {

bool checkHasAllFieldFiles(const std::filesystem::path& snapshot_dir)
{
    const std::vector<std::string> names = {"Ex.gf","Ey.gf","Ez.gf","Hx.gf","Hy.gf","Hz.gf","time.txt"};
    for (const auto& n : names) {
        if (!std::filesystem::exists(snapshot_dir / n)) { return false; }
    }
    return true;
}

double readTimeFile(const std::filesystem::path& snapshot_dir)
{
    const std::filesystem::path time_path = snapshot_dir / "time.txt";
    std::ifstream in(time_path);
    if (!in) { throw std::runtime_error("Missing time.txt in " + snapshot_dir.string()); }
    double t = 0.0;
    in >> t;
    if (!in) { throw std::runtime_error("Invalid time.txt in " + snapshot_dir.string()); }
    return t;
}

std::unique_ptr<GridFunction> loadGridFunction(const std::filesystem::path& file, Mesh& mesh)
{
    std::ifstream in(file, std::ios::binary);
    if (!in) { throw std::runtime_error("Cannot open " + file.string()); }
    auto res = std::make_unique<GridFunction>(&mesh, in);
    return res;
}

TimeToFields buildTimeToFields(const std::string& rank_root, Mesh& mesh)
{
    TimeToFields t2f;
    std::vector<std::pair<double, std::filesystem::path>> time_dirs;
    std::filesystem::path base_dir(rank_root);

    if (!std::filesystem::exists(base_dir) || !std::filesystem::is_directory(base_dir)) {
        throw std::runtime_error("Base directory not found: " + rank_root);
    }

    for (const auto& entry : std::filesystem::directory_iterator(base_dir)) {
        if (!entry.is_directory()) { continue; }
        const std::filesystem::path snap = entry.path();
        if (!checkHasAllFieldFiles(snap)) { continue; }
        const auto time = readTimeFile(snap);
        time_dirs.emplace_back(time, snap);
    }

    std::sort(time_dirs.begin(), time_dirs.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    for (const auto& [t, dir] : time_dirs) {
        TimeFields tfs;
        tfs.Ex = loadGridFunction(dir / "Ex.gf", mesh);
        tfs.Ey = loadGridFunction(dir / "Ey.gf", mesh);
        tfs.Ez = loadGridFunction(dir / "Ez.gf", mesh);
        tfs.Hx = loadGridFunction(dir / "Hx.gf", mesh);
        tfs.Hy = loadGridFunction(dir / "Hy.gf", mesh);
        tfs.Hz = loadGridFunction(dir / "Hz.gf", mesh);
        t2f.emplace(t, std::move(tfs));
    }
    return t2f;
}

}