#pragma once

#include "mfem.hpp"
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <filesystem>

namespace maxwell {

using namespace mfem;

struct TimeFields {
    std::unique_ptr<GridFunction> Ex, Ey, Ez;
    std::unique_ptr<GridFunction> Hx, Hy, Hz;
};

using TimeToFields = std::map<double, TimeFields>;

bool checkHasAllFieldFiles(const std::filesystem::path& snapshot_dir);
double readTimeFile(const std::filesystem::path& snapshot_dir);
std::unique_ptr<GridFunction> loadGridFunction(const std::filesystem::path& file, Mesh& mesh);
TimeToFields buildTimeToFields(const std::string& rank_root, Mesh& mesh);

}