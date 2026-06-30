#pragma once

#include <optional>
#include <string>
#include <vector>

#include "components/FarField.h"

namespace maxwell::driver {

struct MORRCSConfig {
    std::string case_path;
    std::string mesh_path;
    std::string x_dir;
    std::string output_root{"rcs_output"};
    std::string probe_name{"mor_rcs"};
    std::vector<int> tags;
    std::vector<Frequency> frequencies;
    std::vector<SphericalAngles> angles;
    std::optional<double> max_time;
};

void runMORRCSPostProcessing(const MORRCSConfig& cfg);

} // namespace maxwell::driver
