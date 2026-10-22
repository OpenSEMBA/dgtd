#include "RCSSurfacePostProcessor.h"

#include "components/Spherical.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <utility>

#include <omp.h>
#include <unistd.h>

#include "components/RCSSurfaceDftCuda.h"
#include "driver/driver.h"
#include "math/PhysicalConstants.h"
#include "mfemExtension/LinearIntegrators.h"

namespace maxwell {

using namespace mfem;

namespace {

constexpr double kGiB = 1024.0 * 1024.0 * 1024.0;
constexpr double kSystemMemSafetyFraction = 0.5;
constexpr std::uint64_t kMeshHeadroomBytes = 256ull * 1024ull * 1024ull;

static std::vector<std::string> findRankDirs(const std::string& basePath)
{
    std::vector<std::string> paths;
    std::regex pat(R"(rank\d+)");
    for (const auto& e : std::filesystem::directory_iterator(basePath)) {
        if (e.is_directory() &&
            std::regex_match(e.path().filename().string(), pat)) {
            paths.push_back(e.path().string());
        }
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

static void removeDir(const std::string& p)
{
    if (std::filesystem::exists(p)) std::filesystem::remove_all(p);
}

template <typename T>
static T dot3(const std::array<double,3>& a, const std::array<T,3>& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static std::uint64_t readMemAvailableBytes()
{
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo) {
        std::string key;
        std::uint64_t kib = 0;
        std::string unit;
        while (meminfo >> key >> kib >> unit) {
            if (key == "MemAvailable:") {
                return kib * 1024ull;
            }
        }
    }
    const long pages = sysconf(_SC_AVPHYS_PAGES);
    const long pageSize = sysconf(_SC_PAGESIZE);
    if (pages > 0 && pageSize > 0) {
        return static_cast<std::uint64_t>(pages) *
               static_cast<std::uint64_t>(pageSize);
    }
    throw std::runtime_error("Cannot determine available system memory.");
}

static void readGeometryFromStream(std::ifstream& f, RCSSurfaceGeometry& geo)
{
    int32_t hdr[5];
    f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    if (!f) throw std::runtime_error("Failed to read surface_data.bin header.");

    geo.spaceDimension = hdr[0];
    geo.numDofs        = hdr[1];
    geo.numBdrElements = hdr[2];
    geo.numQuadPoints  = hdr[3];
    geo.basisType      = hdr[4];

    const int nqp = geo.numQuadPoints;
    const int sdim = geo.spaceDimension;
    geo.positions.resize(static_cast<std::size_t>(nqp) * sdim);
    geo.normals.resize(static_cast<std::size_t>(nqp) * 3);
    geo.weights.resize(static_cast<std::size_t>(nqp));

    f.read(reinterpret_cast<char*>(geo.positions.data()),
           static_cast<std::streamsize>(geo.positions.size() * sizeof(double)));
    f.read(reinterpret_cast<char*>(geo.normals.data()),
           static_cast<std::streamsize>(geo.normals.size() * sizeof(double)));
    f.read(reinterpret_cast<char*>(geo.weights.data()),
           static_cast<std::streamsize>(geo.weights.size() * sizeof(double)));
    if (!f) throw std::runtime_error("Failed to read surface_data.bin geometry.");
}

static std::uint64_t geometryBytes(const RCSSurfaceGeometry& geo)
{
    return (geo.positions.size() + geo.normals.size() + geo.weights.size())
           * sizeof(double);
}

static std::int64_t snapPayloadBytes(int nDofs)
{
    return static_cast<std::int64_t>(sizeof(double))
           + 6ll * static_cast<std::int64_t>(nDofs) * static_cast<std::int64_t>(sizeof(double));
}

static double elapsedSeconds(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}

static void printElapsedSeconds(const char* label, double seconds)
{
    std::cout << label << std::fixed << std::setprecision(2)
              << seconds << " s\n" << std::defaultfloat;
}

} // namespace

// ------------------------------------------------------------------
// Memory estimate / gate
// ------------------------------------------------------------------

bool RCSSurfacePostProcessor::peekSurfaceFile(
    const std::string& rankPath,
    int everyNSteps,
    RCSSurfaceGeometry& geo,
    std::int64_t& nSnapTotal,
    std::int64_t& nSnapKeptEstimate)
{
    const std::string path = rankPath + "/surface_data.bin";
    const auto fileSize = static_cast<std::uint64_t>(std::filesystem::file_size(path));

    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    readGeometryFromStream(f, geo);

    const auto headerAndGeo =
        static_cast<std::uint64_t>(f.tellg());
    if (fileSize < headerAndGeo) return false;

    const std::int64_t payload = snapPayloadBytes(geo.numDofs);
    if (payload <= 0) return false;

    nSnapTotal = static_cast<std::int64_t>((fileSize - headerAndGeo) / static_cast<std::uint64_t>(payload));
    const int stride = std::max(1, everyNSteps);
    nSnapKeptEstimate = (nSnapTotal + stride - 1) / stride;
    return true;
}

std::uint64_t RCSSurfacePostProcessor::estimateLoadAllBytes(
    std::int64_t nKept, int nDofs, int nFreq)
{
    const std::uint64_t hist =
        static_cast<std::uint64_t>(nKept) * 6ull *
        static_cast<std::uint64_t>(nDofs) * sizeof(double);
    const std::uint64_t ff =
        6ull * static_cast<std::uint64_t>(nFreq) *
        static_cast<std::uint64_t>(nDofs) * sizeof(std::complex<double>);
    return hist + ff + kMeshHeadroomBytes;
}

std::uint64_t RCSSurfacePostProcessor::memoryBudgetBytes(
    const std::optional<double>& ramGateGiB)
{
    if (ramGateGiB.has_value()) {
        if (!(ramGateGiB.value() > 0.0)) {
            throw std::runtime_error("ram_gate must be > 0 (GiB).");
        }
        return static_cast<std::uint64_t>(ramGateGiB.value() * kGiB);
    }
    const auto avail = readMemAvailableBytes();
    return static_cast<std::uint64_t>(kSystemMemSafetyFraction * static_cast<double>(avail));
}

// ------------------------------------------------------------------
// I/O
// ------------------------------------------------------------------

RCSSurfacePostProcessor::RankData
RCSSurfacePostProcessor::readRankData(const std::string& rankPath) const
{
    RankData rd;
    std::ifstream f(rankPath + "/surface_data.bin", std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + rankPath + "/surface_data.bin");

    readGeometryFromStream(f, rd.geometry);

    const int ndofs = rd.geometry.numDofs;
    const int stride = std::max(1, everyNSteps_);
    const std::size_t field_bytes =
        static_cast<std::size_t>(ndofs) * sizeof(double);
    long long snap_index = 0;
    long long kept = 0;
    long long skipped = 0;
    while (f.peek() != EOF) {
        SurfaceSnapshot snap;
        f.read(reinterpret_cast<char*>(&snap.time), sizeof(double));
        if (!f) break;

        const bool keep = (snap_index % stride) == 0;
        if (keep) {
            for (auto* vec : {&snap.Ex, &snap.Ey, &snap.Ez,
                              &snap.Hx, &snap.Hy, &snap.Hz}) {
                vec->resize(ndofs);
                f.read(reinterpret_cast<char*>(vec->data()),
                       static_cast<std::streamsize>(field_bytes));
            }
            if (!f) break;
            rd.snapshots.push_back(std::move(snap));
            ++kept;
        } else {
            // Skip field payload without allocating (critical for large dumps).
            f.seekg(static_cast<std::streamoff>(6 * field_bytes), std::ios::cur);
            if (!f) break;
            ++skipped;
        }
        ++snap_index;
    }
    if (stride > 1) {
        std::cout << "    every_n_steps=" << stride
                  << ": kept " << kept << ", skipped " << skipped
                  << " snapshots while reading\n";
    }
    return rd;
}

RCSSurfacePostProcessor::FreqFields
RCSSurfacePostProcessor::dftFromSnapshots(
    const std::vector<SurfaceSnapshot>& snapshots,
    const std::vector<double>& times,
    const std::vector<double>& normFreqs,
    int nDofs) const
{
    const int nFreq = static_cast<int>(normFreqs.size());
    const int nSnap = static_cast<int>(snapshots.size());
    FreqFields ff(6, std::vector<std::vector<std::complex<double>>>(
        nFreq, std::vector<std::complex<double>>(nDofs, {0, 0})));

    #pragma omp parallel
    for (int t = 0; t < nSnap; ++t) {
        const auto& s = snapshots[t];
        const std::vector<double>* comps[6] = {
            &s.Ex, &s.Ey, &s.Ez, &s.Hx, &s.Hy, &s.Hz };

        #pragma omp for schedule(static) nowait
        for (int fi = 0; fi < nFreq; ++fi) {
            double arg = 2.0 * M_PI * normFreqs[fi] * times[t];
            auto w = std::complex<double>(std::cos(arg), -std::sin(arg));
            for (int c = 0; c < 6; ++c)
                for (int v = 0; v < nDofs; ++v)
                    ff[c][fi][v] += (*comps[c])[v] * w;
        }
    }
    {
        double invN = 1.0 / static_cast<double>(nSnap);
        for (auto& comp : ff)
            for (auto& fv : comp)
                for (auto& v : fv)
                    v *= invN;
    }
    return ff;
}

RCSSurfacePostProcessor::FreqFields
RCSSurfacePostProcessor::dftStreamFromFile(
    const std::string& rankPath,
    const RCSSurfaceGeometry& geo,
    const std::vector<double>& normFreqs,
    std::vector<double>& timesOut,
    int& nSnapOut) const
{
    std::ifstream f(rankPath + "/surface_data.bin", std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + rankPath + "/surface_data.bin");

    RCSSurfaceGeometry discardGeo;
    readGeometryFromStream(f, discardGeo);

    const int nDofs = geo.numDofs;
    const int nFreq = static_cast<int>(normFreqs.size());
    const int stride = std::max(1, everyNSteps_);
    const std::size_t field_bytes =
        static_cast<std::size_t>(nDofs) * sizeof(double);

    FreqFields ff(6, std::vector<std::vector<std::complex<double>>>(
        nFreq, std::vector<std::complex<double>>(nDofs, {0, 0})));

    timesOut.clear();
    long long snap_index = 0;
    long long skipped = 0;

    SurfaceSnapshot snap;
    snap.Ex.resize(nDofs);
    snap.Ey.resize(nDofs);
    snap.Ez.resize(nDofs);
    snap.Hx.resize(nDofs);
    snap.Hy.resize(nDofs);
    snap.Hz.resize(nDofs);

    while (f.peek() != EOF) {
        f.read(reinterpret_cast<char*>(&snap.time), sizeof(double));
        if (!f) break;

        if (maxTime_.has_value() && snap.time > maxTime_.value()) {
            // Times are monotonic in exporter dumps; stop after max_time.
            break;
        }

        const bool keep = (snap_index % stride) == 0;
        if (!keep) {
            f.seekg(static_cast<std::streamoff>(6 * field_bytes), std::ios::cur);
            if (!f) break;
            ++skipped;
            ++snap_index;
            continue;
        }

        for (auto* vec : {&snap.Ex, &snap.Ey, &snap.Ez,
                          &snap.Hx, &snap.Hy, &snap.Hz}) {
            f.read(reinterpret_cast<char*>(vec->data()),
                   static_cast<std::streamsize>(field_bytes));
        }
        if (!f) break;

        const double tNorm = snap.time / physicalConstants::speedOfLight;
        timesOut.push_back(tNorm);

        const std::vector<double>* comps[6] = {
            &snap.Ex, &snap.Ey, &snap.Ez, &snap.Hx, &snap.Hy, &snap.Hz };

        #pragma omp parallel for schedule(static)
        for (int fi = 0; fi < nFreq; ++fi) {
            double arg = 2.0 * M_PI * normFreqs[fi] * tNorm;
            auto w = std::complex<double>(std::cos(arg), -std::sin(arg));
            for (int c = 0; c < 6; ++c)
                for (int v = 0; v < nDofs; ++v)
                    ff[c][fi][v] += (*comps[c])[v] * w;
        }

        ++snap_index;
    }

    nSnapOut = static_cast<int>(timesOut.size());
    if (nSnapOut == 0) {
        throw std::runtime_error("No snapshots found while streaming surface_data.bin.");
    }

    {
        double invN = 1.0 / static_cast<double>(nSnapOut);
        for (auto& comp : ff)
            for (auto& fv : comp)
                for (auto& v : fv)
                    v *= invN;
    }

    if (stride > 1) {
        std::cout << "    every_n_steps=" << stride
                  << ": kept " << nSnapOut << ", skipped " << skipped
                  << " snapshots while streaming\n";
    }
    return ff;
}

bool RCSSurfacePostProcessor::tryCudaStreamDft(
    const std::string& rankPath,
    int nDofs,
    const std::vector<double>& normFreqs,
    std::vector<double>& timesOut,
    int& nSnapOut,
    FreqFields& ffOut) const
{
#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (!mfem::Device::Allows(mfem::Backend::CUDA)) {
        return false;
    }
    return rcsCudaDftStreamFile(
        rankPath + "/surface_data.bin",
        nDofs,
        everyNSteps_,
        maxTime_,
        normFreqs,
        timesOut,
        nSnapOut,
        ffOut);
#else
    (void)rankPath;
    (void)nDofs;
    (void)normFreqs;
    (void)timesOut;
    (void)nSnapOut;
    (void)ffOut;
    return false;
#endif
}

// ------------------------------------------------------------------
// Plane-wave helpers
// ------------------------------------------------------------------

PlaneWaveData RCSSurfacePostProcessor::extractPlaneWaveData(
    const std::string& jsonPath) const
{
    std::string meshDir = std::filesystem::path(jsonPath).parent_path().string();
    return buildPlaneWaveData(driver::parseJSONfile(jsonPath), meshDir);
}

std::vector<double> RCSSurfacePostProcessor::computeIncidentPowerSpectrum(
    const PlaneWaveData& pw,
    const std::vector<double>& times,
    const std::vector<Frequency>& freqs) const
{
    auto envelope = evaluateGaussianVector(
        const_cast<std::vector<double>&>(times), pw.spread, pw.mean);

    // Apply carrier modulation when using a modulated Gaussian.
    if (pw.isModulated()) {
        for (size_t t = 0; t < times.size(); ++t) {
            double carrier_arg = 2.0 * M_PI * pw.frequency * (times[t] - std::abs(pw.mean));
            envelope[t] *= std::cos(carrier_arg);
        }
    }

    std::vector<double> power(freqs.size(), 0.0);
    for (size_t fi = 0; fi < freqs.size(); ++fi) {
        std::complex<double> val(0.0, 0.0);
        for (size_t t = 0; t < times.size(); ++t) {
            double arg = 2.0 * M_PI * freqs[fi] * times[t];
            val += envelope[t] * std::complex<double>(std::cos(arg), -std::sin(arg));
        }
        val /= static_cast<double>(times.size());
        power[fi] = std::norm(val) / (2.0 * physicalConstants::freeSpaceImpedance);
    }
    return power;
}

// ------------------------------------------------------------------
// Determine FEC order by matching ndofs
// ------------------------------------------------------------------

static int determineFECOrder(ParMesh& pmesh, int nDofs)
{
    int meshDim = pmesh.Dimension();
    for (int p = 0; p <= 10; ++p) {
        DG_FECollection testFec(p, meshDim);
        ParFiniteElementSpace testFes(&pmesh, &testFec);
        if (testFes.GetNDofs() == nDofs) return p;
    }
    throw std::runtime_error("Could not determine FEC order from numDofs.");
}

// ------------------------------------------------------------------
// Core computation
// ------------------------------------------------------------------

void RCSSurfacePostProcessor::computeAndWriteResults(
    const std::string& dataPath,
    const std::string& jsonPath,
    std::vector<Frequency>& frequencies,
    const std::vector<SphericalAngles>& angles)
{
    auto rankPaths = findRankDirs(dataPath);
    if (rankPaths.empty())
        throw std::runtime_error("No rank folders in " + dataPath);

    std::cout << "[RCS] Post-processing started.\n"
              << "  Data path  : " << dataPath << "\n"
              << "  Ranks      : " << rankPaths.size()
              << "   Frequencies: " << frequencies.size()
              << "   Angles: " << angles.size() << "\n";

    // Rescale frequencies from Hz to normalised (f / c_SI).
    std::vector<double> normFreqs(frequencies.size());
    for (size_t i = 0; i < frequencies.size(); ++i)
        normFreqs[i] = frequencies[i] / physicalConstants::speedOfLight_SI;

    const int nFreq = static_cast<int>(normFreqs.size());
    const std::uint64_t budget = memoryBudgetBytes(ramGateGiB_);
    if (ramGateGiB_.has_value()) {
        std::cout << "  ram_gate   : " << ramGateGiB_.value() << " GiB ("
                  << std::fixed << std::setprecision(2)
                  << (budget / kGiB) << " GiB hard budget)\n"
                  << std::defaultfloat;
    } else {
        std::cout << "  RAM budget : " << std::fixed << std::setprecision(2)
                  << (budget / kGiB) << " GiB "
                  << "(50% of MemAvailable)\n"
                  << std::defaultfloat;
    }

    // Initialise output maps.
    for (const auto& ang : angles)
        for (const auto& f : normFreqs) {
            farFieldData_[ang][f] = 0.0;
            rcsData_[ang][f] = 0.0;
        }

    // Accumulators for coherent summation of complex amplitudes across ranks.
    // Key: (angle, frequency), Value: {N_theta, N_phi, L_theta, L_phi}
    std::map<std::pair<SphericalAngles, double>, std::array<std::complex<double>, 4>> coherentSum;
    for (const auto& ang : angles)
        for (const auto& f : normFreqs)
            coherentSum[{ang, f}] = {0, 0, 0, 0};

    PlaneWaveData pw(0.0, 0.0);
    std::vector<double> incidentPower;
    int spaceDim = 0;
    bool firstRank = true;
    size_t rankIdx = 0;

    for (const auto& rp : rankPaths) {
        ++rankIdx;

        RCSSurfaceGeometry peekGeo;
        std::int64_t nSnapTotal = 0;
        std::int64_t nSnapKeptEst = 0;
        if (!peekSurfaceFile(rp, everyNSteps_, peekGeo, nSnapTotal, nSnapKeptEst)) {
            throw std::runtime_error("Cannot peek surface_data.bin in " + rp);
        }

        const int nDofs = peekGeo.numDofs;
        spaceDim = peekGeo.spaceDimension;
        const std::uint64_t estBytes =
            estimateLoadAllBytes(nSnapKeptEst, nDofs, nFreq) + geometryBytes(peekGeo);
        const bool useStream = estBytes > budget;
#ifdef SEMBA_DGTD_ENABLE_CUDA
        const bool cudaCandidate = mfem::Device::Allows(mfem::Backend::CUDA);
#else
        const bool cudaCandidate = false;
#endif
        const char* readMode = cudaCandidate
            ? "CUDA streaming"
            : (useStream ? "Streaming" : "Load-all");

        std::cout << "\n  [Rank " << rankIdx << "/" << rankPaths.size()
                  << "] " << readMode
                  << " surface data..."
                  << " (est " << std::fixed << std::setprecision(2)
                  << (estBytes / kGiB) << " GiB vs budget "
                  << (budget / kGiB) << " GiB; "
                  << nSnapKeptEst << " snaps est, " << nDofs << " DOFs)\n"
                  << std::defaultfloat;

        if (spaceDim == 2)
            for (const auto& a : angles)
                if (std::abs(a.theta - M_PI_2) > 1e-8)
                    throw std::runtime_error("2D RCS requires theta = pi/2.");

        FreqFields ff;
        std::vector<double> times;
        int nSnap = 0;
        int basisType = peekGeo.basisType;

        const auto dftT0 = std::chrono::steady_clock::now();
        const bool usedCuda = tryCudaStreamDft(rp, nDofs, normFreqs, times, nSnap, ff);
        if (usedCuda) {
            std::cout << "    Stream DFT : " << nSnap << " snapshots x "
                      << nFreq << " frequencies x " << nDofs << " DOFs done.\n";
            printElapsedSeconds("    DFT time   : ", elapsedSeconds(dftT0));
            if (maxTime_.has_value()) {
                std::cout << "    Time filter: applied during stream (maxTime = "
                          << std::fixed << std::setprecision(5)
                          << maxTime_.value() << ").\n" << std::defaultfloat;
            }
        } else if (useStream) {
            ff = dftStreamFromFile(rp, peekGeo, normFreqs, times, nSnap);
            std::cout << "    Stream DFT : " << nSnap << " snapshots x "
                      << nFreq << " frequencies x " << nDofs << " DOFs done.\n";
            printElapsedSeconds("    DFT time   : ", elapsedSeconds(dftT0));
            if (maxTime_.has_value()) {
                std::cout << "    Time filter: applied during stream (maxTime = "
                          << std::fixed << std::setprecision(5)
                          << maxTime_.value() << ").\n" << std::defaultfloat;
            }
        } else {
            std::cout << "    Reading..." << std::flush;
            auto rd = readRankData(rp);
            basisType = rd.geometry.basisType;
            // Move — do not deep-copy the snapshot history (peak RAM).
            std::vector<SurfaceSnapshot> snapshots = std::move(rd.snapshots);
            std::cout << " done. (" << spaceDim << "D, " << nDofs
                      << " DOFs, " << snapshots.size() << " snapshots)\n";

            if (maxTime_.has_value()) {
                const size_t before = snapshots.size();
                snapshots.erase(
                    std::remove_if(snapshots.begin(), snapshots.end(),
                        [this](const SurfaceSnapshot& s) {
                            return s.time > maxTime_.value();
                        }),
                    snapshots.end());
                if (snapshots.empty()) {
                    throw std::runtime_error(
                        "No snapshots found within the specified maxTime.");
                }
                std::cout << "    Time filter: kept " << snapshots.size() << "/"
                          << before << " snapshots (maxTime = " << std::fixed
                          << std::setprecision(5) << maxTime_.value() << ").\n"
                          << std::defaultfloat;
            }
            nSnap = static_cast<int>(snapshots.size());
            times.resize(static_cast<std::size_t>(nSnap));
            for (int i = 0; i < nSnap; ++i)
                times[i] = snapshots[i].time / physicalConstants::speedOfLight;

            std::cout << "    DFT        : " << nSnap << " snapshots x "
                      << nFreq << " frequencies x " << nDofs << " DOFs..."
                      << std::flush;
            ff = dftFromSnapshots(snapshots, times, normFreqs, nDofs);
            snapshots.clear();
            snapshots.shrink_to_fit();
            std::cout << " done.\n";
            printElapsedSeconds("    DFT time   : ", elapsedSeconds(dftT0));
        }

        if (firstRank) {
            pw = extractPlaneWaveData(jsonPath);
            std::cout << "    Plane-wave : spread=" << pw.spread
                      << ", mean=" << pw.mean;
            if (pw.isModulated())
                std::cout << " [modulated at f=" << std::scientific
                          << std::setprecision(2)
                          << pw.frequency * physicalConstants::speedOfLight_SI << " Hz]";
            std::cout << ".\n" << std::defaultfloat;
            incidentPower = computeIncidentPowerSpectrum(pw, times, normFreqs);
            std::cout << "    Incident power spectrum computed for "
                      << normFreqs.size() << " frequencies.\n";
            firstRank = false;
        }

        // --- Load mesh and build FES for this rank ---
        std::cout << "    Loading mesh..." << std::flush;
        auto mesh = Mesh::LoadFromFile(rp + "/mesh", 1, 0);
        auto pmesh = ParMesh(MPI_COMM_WORLD, mesh);
        int order = determineFECOrder(pmesh, nDofs);
        DG_FECollection fec(order, pmesh.Dimension(), basisType);
        ParFiniteElementSpace fes(&pmesh, &fec);
        std::cout << " done. (FEC order: " << order << ")\n";

        // --- For each frequency and angle, compute far-field potentials ---
        const auto farT0 = std::chrono::steady_clock::now();
        std::cout << "    Far-field integration: " << nFreq
                  << " frequencies x " << angles.size() << " angles..." << std::flush;
        for (int fi = 0; fi < nFreq; ++fi) {
            const double freq = normFreqs[fi];

            for (const auto& ang : angles) {
                // Build phase-term function coefficients.
                std::unique_ptr<FunctionCoefficient> fcR, fcI;
                fcR = buildFC(spaceDim, freq, ang, true);
                fcI = buildFC(spaceDim, freq, ang, false);

                // For each spatial direction d, assemble linear forms that
                // compute:  lf[i] = integral n[d] * fc * shape_i dS
                // over NTF boundary faces.
                std::array<std::complex<double>, 3> N_vec = {0, 0, 0};
                std::array<std::complex<double>, 3> L_vec = {0, 0, 0};

                for (int dir = X; dir <= Z; ++dir) {
                    auto lfR = assembleLinearForm(*fcR, fes, dir);
                    auto lfI = assembleLinearForm(*fcI, fes, dir);

                    // For each field component c, compute complex integral:
                    //   I_{dir,c} = integral n[dir] * fc_complex * Field_c dS
                    // where fc_complex = fcR + j*fcI
                    // I_{dir,c} = sum_i (lfR[i] + j*lfI[i]) * Field_c_dof[i]
                    // This equals: sum_i lfR[i]*Re(F_i) - lfI[i]*Im(F_i)
                    //            + j*(lfR[i]*Im(F_i) + lfI[i]*Re(F_i))

                    std::array<std::complex<double>, 3> Hint, Eint;
                    for (int c = 0; c < 3; ++c) {
                        auto& Hf = ff[3 + c][fi];
                        auto& Ef = ff[c][fi];
                        double rH = 0, iH = 0, rE = 0, iE = 0;
                        for (int v = 0; v < nDofs; ++v) {
                            rH += lfR->Elem(v) * Hf[v].real() - lfI->Elem(v) * Hf[v].imag();
                            iH += lfR->Elem(v) * Hf[v].imag() + lfI->Elem(v) * Hf[v].real();
                            rE += lfR->Elem(v) * Ef[v].real() - lfI->Elem(v) * Ef[v].imag();
                            iE += lfR->Elem(v) * Ef[v].imag() + lfI->Elem(v) * Ef[v].real();
                        }
                        Hint[c] = {rH, iH};
                        Eint[c] = {rE, iE};
                    }

                    // Accumulate cross-product terms for J = n x H, M = -n x E.
                    // (n x H)_i = eps_{ijk} n_j H_k  =>  for dir=j:
                    int j = dir;
                    int i1 = (j + 1) % 3, k1 = (j + 2) % 3;
                    int i2 = (j + 2) % 3, k2 = (j + 1) % 3;
                    N_vec[i1] += Hint[k1];
                    N_vec[i2] -= Hint[k2];
                    L_vec[i1] -= Eint[k1];   // M = -n x E
                    L_vec[i2] += Eint[k2];
                }

                // Project onto spherical components.
                auto th = thetaHat(ang.theta, ang.phi);
                auto ph = phiHat(ang.phi);

                auto N_theta = dot3(th, N_vec);
                auto N_phi   = dot3(ph, N_vec);
                auto L_theta = dot3(th, L_vec);
                auto L_phi   = dot3(ph, L_vec);

                // Accumulate complex amplitudes coherently across ranks.
                auto& acc = coherentSum[{ang, freq}];
                acc[0] += N_theta;
                acc[1] += N_phi;
                acc[2] += L_theta;
                acc[3] += L_phi;
            }
        }
        std::cout << " done.\n";
        printElapsedSeconds("    Far-field time: ", elapsedSeconds(farT0));
    }

    // After all ranks processed, compute far-field power from coherently-summed amplitudes.
    std::cout << "\n[RCS] Computing far-field radiation potentials..." << std::flush;
    for (int fi = 0; fi < nFreq; ++fi) {
        const double freq = normFreqs[fi];
        const double k = 2.0 * M_PI * freq;
        double Z0 = physicalConstants::freeSpaceImpedance;

        for (const auto& ang : angles) {
            const auto& acc = coherentSum[{ang, freq}];
            auto N_theta = acc[0];
            auto N_phi   = acc[1];
            auto L_theta = acc[2];
            auto L_phi   = acc[3];

            double potRad;
            if (spaceDim == 2) {
                // 2D NTFF formula derived from the scalar Green's function for E_y (TE):
                //   E_y^scat ~ k/(4j) * sqrt(2/(pi*k*r)) * e^{-j(kr-pi/4)} * integral J_y e^{jk r_hat.r'} dl
                // sigma_2D = 2*pi*r * |E_scat|^2/|E_inc|^2 = k/4 * |N_phi|^2 / |E_inc|^2
                // With potRad = k^2/(32*Z0) * |N|^2 and ct = 4/k / P_inc:
                //   sigma = (4/k) * k^2/(32*Z0) * 2*Z0 / |E_inc|^2 = k/4 * |N|^2 / |E_inc|^2  ✓
                potRad = k * k / (32.0 * Z0) *
                    (std::norm(L_phi + Z0 * N_theta) +
                     std::norm(L_theta - Z0 * N_phi));
            } else {
                potRad = k * k / (32.0 * M_PI * M_PI * Z0) *
                    (std::norm(L_phi + Z0 * N_theta) +
                     std::norm(L_theta - Z0 * N_phi));
            }
            farFieldData_[ang][freq] = potRad;
        }
    }

    std::cout << " done.\n";

    // Compute RCS from accumulated far-field data.
    std::cout << "[RCS] Computing RCS..." << std::flush;
    for (int fi = 0; fi < nFreq; ++fi) {
        double k = 2.0 * M_PI * normFreqs[fi];
        for (const auto& ang : angles) {
            double ct;
            if (spaceDim == 2)
                ct = 4.0 / (k * incidentPower[fi]);
            else
                ct = 4.0 * M_PI / incidentPower[fi];
            rcsData_[ang][normFreqs[fi]] = ct * farFieldData_[ang][normFreqs[fi]];
        }
    }

    std::cout << " done.\n";

    // Write merged output files.
    std::cout << "[RCS] Writing results to '" << dataPath << "'..." << std::flush;
    removeDir(dataPath + "/farfield");
    removeDir(dataPath + "/rcs");
    std::filesystem::create_directories(dataPath + "/farfield");
    std::filesystem::create_directories(dataPath + "/rcs");

    for (const auto& ang : angles) {
        std::string sfx = "Th_" + std::to_string(ang.theta) +
                           "_Phi_" + std::to_string(ang.phi) + "_dgtd.dat";
        {
            std::ofstream out(dataPath + "/farfield/farfieldData_" + sfx);
            out << "Theta (rad) // Phi (rad) // Frequency (Hz) // pot // normalization_term\n";
            for (const auto& f : normFreqs) {
                double lam = physicalConstants::speedOfLight / f;
                double norm = (spaceDim == 2) ? lam : lam * lam;
                out << ang.theta << " " << ang.phi << " "
                    << f * physicalConstants::speedOfLight_SI << " "
                    << farFieldData_[ang][f] << " " << norm << "\n";
            }
        }
        {
            std::ofstream out(dataPath + "/rcs/rcsData_" + sfx);
            out << "Theta (rad) // Phi (rad) // Frequency (Hz) // rcs // normalization_term\n";
            for (const auto& f : normFreqs) {
                double lam = physicalConstants::speedOfLight / f;
                double norm = (spaceDim == 2) ? lam : lam * lam;
                out << ang.theta << " " << ang.phi << " "
                    << f * physicalConstants::speedOfLight_SI << " "
                    << rcsData_[ang][f] << " " << norm << "\n";
            }
        }
    }
    std::cout << " done. (" << angles.size() << " angle(s)).\n";
    std::cout << "[RCS] Post-processing complete.\n";
}

// ------------------------------------------------------------------
// Constructor
// ------------------------------------------------------------------

RCSSurfacePostProcessor::RCSSurfacePostProcessor(
    const std::string& dataPath,
    const std::string& jsonPath,
    std::vector<Frequency>& frequencies,
    const std::vector<SphericalAngles>& angles,
    const std::optional<double>& maxTime,
    int everyNSteps,
    const std::optional<double>& ramGateGiB)
    : maxTime_(maxTime),
      everyNSteps_(std::max(1, everyNSteps)),
      ramGateGiB_(ramGateGiB)
{
    computeAndWriteResults(dataPath, jsonPath, frequencies, angles);
}

} // namespace maxwell
