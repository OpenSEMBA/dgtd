#pragma once

#include <complex>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "components/FarField.h"

namespace maxwell {

struct RCSSurfaceGeometry {
    int spaceDimension;
    int numDofs;
    int numBdrElements;
    int numQuadPoints;
    int basisType;

    std::vector<double> positions;  // [numQuadPoints * spaceDimension]
    std::vector<double> normals;    // [numQuadPoints * 3]
    std::vector<double> weights;    // [numQuadPoints]
};

struct SurfaceSnapshot {
    double time;
    std::vector<double> Ex, Ey, Ez;
    std::vector<double> Hx, Hy, Hz;
};

/**
 * Reads surface data exported by RCSSurfaceExporter and computes RCS.
 *
 * The algorithm:
 *   1. Read geometry (positions, normals, weights) and time-domain snapshots.
 *   2. For each requested frequency, DFT the fields.
 *   3. Compute equivalent surface currents via n x H and -n x E.
 *   4. Compute far-field radiation potentials via surface integration with
 *      the phase kernel exp(+j k r'.r_hat).
 *   5. Normalise to obtain RCS.
 *
 * Memory: estimates peak RAM for load-all vs a budget (optional `ram_gate` GiB,
 * else ~50% of MemAvailable). If the dump fits, load-all + move; otherwise
 * stream time samples into the frequency-domain accumulator without retaining
 * the full history (same Nt / times ⇒ near-equivalent DFT).
 *
 * Works identically for 2D and 3D; the only difference is the RCS
 * normalisation constant and the 2D restriction theta=pi/2.
 */
class RCSSurfacePostProcessor {
public:
    RCSSurfacePostProcessor(
        const std::string& dataPath,
        const std::string& jsonPath,
        std::vector<Frequency>& frequencies,
        const std::vector<SphericalAngles>& angles,
        const std::optional<double>& maxTime = std::nullopt,
        int everyNSteps = 1,
        const std::optional<double>& ramGateGiB = std::nullopt);

    double getRCS(const SphericalAngles& angles, const Frequency& freq) const {
        return rcsData_.at(angles).at(freq);
    }

    double getFarField(const SphericalAngles& angles, const Frequency& freq) const {
        return farFieldData_.at(angles).at(freq);
    }

private:
    struct RankData {
        RCSSurfaceGeometry geometry;
        std::vector<SurfaceSnapshot> snapshots;
    };

    using FreqFields = std::vector<std::vector<std::vector<std::complex<double>>>>;

    RankData readRankData(const std::string& rankPath) const;

    /** Peek header + estimate kept snapshot count (stride only; max_time ignored). */
    static bool peekSurfaceFile(
        const std::string& rankPath,
        int everyNSteps,
        RCSSurfaceGeometry& geo,
        std::int64_t& nSnapTotal,
        std::int64_t& nSnapKeptEstimate);

    static std::uint64_t estimateLoadAllBytes(
        std::int64_t nKept, int nDofs, int nFreq);

    static std::uint64_t memoryBudgetBytes(const std::optional<double>& ramGateGiB);

    FreqFields dftFromSnapshots(
        const std::vector<SurfaceSnapshot>& snapshots,
        const std::vector<double>& times,
        const std::vector<double>& normFreqs,
        int nDofs) const;

    /** Stream surface_data.bin into DFT + times; apply stride and max_time. */
    FreqFields dftStreamFromFile(
        const std::string& rankPath,
        const RCSSurfaceGeometry& geo,
        const std::vector<double>& normFreqs,
        std::vector<double>& timesOut,
        int& nSnapOut) const;

    PlaneWaveData extractPlaneWaveData(const std::string& jsonPath) const;

    std::vector<double> computeIncidentPowerSpectrum(
        const PlaneWaveData& pw,
        const std::vector<double>& times,
        const std::vector<Frequency>& frequencies) const;

    void computeAndWriteResults(
        const std::string& dataPath,
        const std::string& jsonPath,
        std::vector<Frequency>& frequencies,
        const std::vector<SphericalAngles>& angles);

    std::map<SphericalAngles, Freq2Value> rcsData_;
    std::optional<double> maxTime_;
    int everyNSteps_{1};  // keep snapshot i if i % everyNSteps_ == 0 (applied while reading)
    std::optional<double> ramGateGiB_;  // optional hard budget (GiB); omit → system MemAvailable
    std::map<SphericalAngles, Freq2Value> farFieldData_;
};

} // namespace maxwell
