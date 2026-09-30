#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <mfem.hpp>

namespace maxwell::driver {

/// Contract with mor_dgtd operator dumps (see docs/mor2rcs.md).
struct MorUrMeta {
    std::int64_t N = -1;   // rows = 6 * ndofs
    std::int64_t r = -1;   // cols = reduced rank
    std::string dtype{"float64"};
    std::string layout{"colmajor"}; // only colmajor float64 supported
    std::string ur_file{"Ur.bin"};
    std::optional<double> dt;
    std::optional<double> final_time;
    std::string mor_case;
};

MorUrMeta loadMorUrMeta(const std::string& meta_path);

/// Load dense Ur (N×r) from a raw binary file. Bytes = N*r*sizeof(double).
/// layout must be column-major (Fortran): entry (i,j) at offset i + j*N.
std::vector<double> loadUrBinary(const std::string& ur_path,
                                 std::int64_t N,
                                 std::int64_t r,
                                 const std::string& layout = "colmajor",
                                 const std::string& dtype = "float64");

/// Collect xr_<k> (or x_<k>) ASCII snapshots sorted by index.
std::vector<std::pair<int, std::string>> collectIndexedSnapshots(
    const std::string& dir,
    const std::string& prefix);

/// ASCII snapshot: time, size, then `size` doubles (same header as MOR x_/xr_).
std::pair<double, mfem::Vector> loadAsciiSnapshot(const std::string& file_path,
                                                  int expected_size);

/// x = Ur * xr with Ur column-major (N×r). OpenMP over rows when available.
void urMatVec(const std::vector<double>& Ur,
              std::int64_t N,
              std::int64_t r,
              const mfem::Vector& xr,
              mfem::Vector& x);

} // namespace maxwell::driver
