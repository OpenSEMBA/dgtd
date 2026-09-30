#pragma once

#include <optional>
#include <string>
#include <vector>

#include "components/FarField.h"

namespace maxwell::driver {

struct MORRCSConfig {
    std::string case_path;
    std::string mesh_path;
    std::string x_dir;      // legacy ASCII full-order snapshots (x_k)
    std::string ur_path;    // optional: dense Ur.bin (N×r, float64 col-major)
    std::string xr_dir;     // optional: reduced snapshots xr_k
    std::string meta_path;  // optional: meta.json next to Ur (auto if empty)
    std::string dump_x_dir; // optional debug: write reconstructed x_k ASCII
    std::string output_root{"rcs_output"};
    std::string probe_name{"mor_rcs"};
    std::vector<int> tags;
    std::vector<Frequency> frequencies;
    std::vector<SphericalAngles> angles;
    std::optional<double> max_time;

    bool useUrReconstruct() const
    {
        return !ur_path.empty() && !xr_dir.empty();
    }
};

void runMORRCSPostProcessing(const MORRCSConfig& cfg);

} // namespace maxwell::driver
