#include "MorUrReconstructor.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <stdexcept>

#include <nlohmann/json.hpp>

#ifdef SEMBA_DGTD_ENABLE_OPENMP
#include <omp.h>
#endif

namespace maxwell::driver {

namespace {

using json = nlohmann::json;

void requireFinite(const mfem::Vector& v, const std::string& what)
{
    for (int i = 0; i < v.Size(); ++i) {
        if (!std::isfinite(v[i])) {
            throw std::runtime_error(what + " contains non-finite value at index " +
                                     std::to_string(i));
        }
    }
}

} // namespace

MorUrMeta loadMorUrMeta(const std::string& meta_path)
{
    std::ifstream f(meta_path);
    if (!f) {
        throw std::runtime_error("Cannot open MOR meta JSON: " + meta_path);
    }
    const json j = json::parse(f);

    MorUrMeta meta;
    meta.N = j.at("N").get<std::int64_t>();
    meta.r = j.at("r").get<std::int64_t>();
    if (j.contains("dtype")) {
        meta.dtype = j.at("dtype").get<std::string>();
    }
    if (j.contains("layout")) {
        meta.layout = j.at("layout").get<std::string>();
    }
    if (j.contains("ur_file")) {
        meta.ur_file = j.at("ur_file").get<std::string>();
    }
    if (j.contains("dt") && !j["dt"].is_null()) {
        meta.dt = j.at("dt").get<double>();
    }
    if (j.contains("final_time") && !j["final_time"].is_null()) {
        meta.final_time = j.at("final_time").get<double>();
    }
    if (j.contains("mor_case")) {
        meta.mor_case = j.at("mor_case").get<std::string>();
    }

    if (meta.N <= 0 || meta.r <= 0) {
        throw std::runtime_error("Invalid N/r in MOR meta: " + meta_path);
    }
    return meta;
}

std::vector<double> loadUrBinary(const std::string& ur_path,
                                 std::int64_t N,
                                 std::int64_t r,
                                 const std::string& layout,
                                 const std::string& dtype)
{
    if (dtype != "float64" && dtype != "double") {
        throw std::runtime_error("Unsupported Ur dtype '" + dtype +
                                 "' (only float64/double).");
    }
    if (layout != "colmajor" && layout != "column-major" && layout != "fortran") {
        throw std::runtime_error("Unsupported Ur layout '" + layout +
                                 "' (only colmajor).");
    }
    if (N <= 0 || r <= 0) {
        throw std::runtime_error("Ur dimensions must be positive.");
    }

    const std::uint64_t expected_bytes =
        static_cast<std::uint64_t>(N) * static_cast<std::uint64_t>(r) * sizeof(double);

    std::ifstream ifs(ur_path, std::ios::binary | std::ios::ate);
    if (!ifs) {
        throw std::runtime_error("Cannot open Ur binary: " + ur_path);
    }
    const auto file_bytes = static_cast<std::uint64_t>(ifs.tellg());
    if (file_bytes != expected_bytes) {
        throw std::runtime_error(
            "Ur.bin size mismatch for " + ur_path + ": got " +
            std::to_string(file_bytes) + " bytes, expected " +
            std::to_string(expected_bytes) + " (= N*r*8 with N=" +
            std::to_string(N) + ", r=" + std::to_string(r) + ").");
    }
    ifs.seekg(0, std::ios::beg);

    std::vector<double> Ur(static_cast<std::size_t>(N * r));
    ifs.read(reinterpret_cast<char*>(Ur.data()),
             static_cast<std::streamsize>(expected_bytes));
    if (!ifs) {
        throw std::runtime_error("Failed reading Ur binary: " + ur_path);
    }
    return Ur;
}

std::vector<std::pair<int, std::string>> collectIndexedSnapshots(
    const std::string& dir,
    const std::string& prefix)
{
    // Accept xr_0 / x_0 and optional .txt (MOR writes xr_<k>.txt).
    const std::regex re("^" + prefix + R"(_(\d+)(?:\.txt)?$)");
    std::vector<std::pair<int, std::string>> files;

    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        std::smatch m;
        const auto name = entry.path().filename().string();
        if (std::regex_match(name, m, re)) {
            files.emplace_back(std::stoi(m[1].str()), entry.path().string());
        }
    }

    std::sort(files.begin(), files.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    if (files.empty()) {
        throw std::runtime_error("No " + prefix + "_k snapshots found in: " + dir);
    }
    return files;
}

std::pair<double, mfem::Vector> loadAsciiSnapshot(const std::string& file_path,
                                                  int expected_size)
{
    std::ifstream ifs(file_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open snapshot: " + file_path);
    }

    double time = 0.0;
    long long n = -1;
    if (!(ifs >> time)) {
        throw std::runtime_error("Invalid time header in: " + file_path);
    }
    if (!(ifs >> n)) {
        throw std::runtime_error("Invalid size header in: " + file_path);
    }
    if (n != expected_size) {
        throw std::runtime_error(
            "Snapshot size mismatch in " + file_path + ". Expected " +
            std::to_string(expected_size) + ", got " + std::to_string(n));
    }

    mfem::Vector state(expected_size);
    for (int i = 0; i < expected_size; ++i) {
        if (!(ifs >> state[i])) {
            throw std::runtime_error(
                "Unexpected EOF while reading state values in: " + file_path);
        }
    }
    return {time, std::move(state)};
}

void urMatVec(const std::vector<double>& Ur,
              std::int64_t N,
              std::int64_t r,
              const mfem::Vector& xr,
              mfem::Vector& x)
{
    if (xr.Size() != static_cast<int>(r)) {
        throw std::runtime_error(
            "xr size " + std::to_string(xr.Size()) + " != Ur cols r=" +
            std::to_string(r));
    }
    if (x.Size() != static_cast<int>(N)) {
        x.SetSize(static_cast<int>(N));
    }
    if (static_cast<std::int64_t>(Ur.size()) != N * r) {
        throw std::runtime_error("Ur buffer size does not match N*r.");
    }

    const double* A = Ur.data();
    const double* z = xr.GetData();
    double* y = x.GetData();

#ifdef SEMBA_DGTD_ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < N; ++i) {
        double sum = 0.0;
        for (std::int64_t j = 0; j < r; ++j) {
            sum += A[i + j * N] * z[j];
        }
        y[i] = sum;
    }

    requireFinite(x, "Reconstructed full-order state x=Ur*xr");
}

} // namespace maxwell::driver
