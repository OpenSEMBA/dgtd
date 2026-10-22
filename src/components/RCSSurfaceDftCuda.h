#pragma once

#include <complex>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace maxwell {

// Host spectra layout matches RCSSurfacePostProcessor::FreqFields:
// ff[component][frequency][dof], component order Ex, Ey, Ez, Hx, Hy, Hz.
using RcsFreqFields =
    std::vector<std::vector<std::vector<std::complex<double>>>>;

struct RcsCudaDftTiles {
    int dofTile{0};
    int freqTile{0};
    bool ok{false};
};

inline int rcsCudaDftPassCount(int n, int tile)
{
    if (tile <= 0 || n <= 0) {
        return 0;
    }
    return (n + tile - 1) / tile;
}

// Bytes for one complex accumulator entry of all six field components.
inline constexpr std::uint64_t rcsCudaDftAccBytesPerFreqDof()
{
    return 6ull * sizeof(std::complex<double>);
}

// Choose resident accumulator tiles that fit in `accBudgetBytes`.
// Prefer one tile covering every frequency and DOF. Otherwise tile DOFs and
// keep every frequency. Split frequencies only when a single DOF across all
// frequencies does not fit. Each extra tile is another pass over the snapshots.
inline RcsCudaDftTiles chooseRcsCudaDftTiles(
    int nFreq, int nDofs, std::uint64_t accBudgetBytes)
{
    RcsCudaDftTiles tiles;
    const std::uint64_t bytesPer = rcsCudaDftAccBytesPerFreqDof();
    if (nFreq <= 0 || nDofs <= 0 || accBudgetBytes < bytesPer) {
        return tiles;
    }

    auto fits = [&](int nFreqTile, int nDofTile) -> bool {
        if (nFreqTile <= 0 || nDofTile <= 0) {
            return false;
        }
        if (static_cast<std::uint64_t>(nFreqTile) > accBudgetBytes / bytesPer) {
            return false;
        }
        const std::uint64_t maxDof =
            accBudgetBytes / (bytesPer * static_cast<std::uint64_t>(nFreqTile));
        return static_cast<std::uint64_t>(nDofTile) <= maxDof;
    };

    if (fits(nFreq, nDofs)) {
        tiles.dofTile = nDofs;
        tiles.freqTile = nFreq;
        tiles.ok = true;
        return tiles;
    }

    if (fits(nFreq, 1)) {
        std::uint64_t maxDof =
            accBudgetBytes / (bytesPer * static_cast<std::uint64_t>(nFreq));
        if (maxDof > static_cast<std::uint64_t>(nDofs)) {
            maxDof = static_cast<std::uint64_t>(nDofs);
        }
        int dofTile = static_cast<int>(maxDof);
        if (dofTile > 256) {
            dofTile -= dofTile % 128;
        }
        if (dofTile < 1) {
            return tiles;
        }
        tiles.dofTile = dofTile;
        tiles.freqTile = nFreq;
        tiles.ok = true;
        return tiles;
    }

    int modest = nDofs < 4096 ? nDofs : 4096;
    while (modest > 1 && !fits(1, modest)) {
        modest /= 2;
    }
    if (!fits(1, modest)) {
        return tiles;
    }
    std::uint64_t maxFreq =
        accBudgetBytes / (bytesPer * static_cast<std::uint64_t>(modest));
    if (maxFreq > static_cast<std::uint64_t>(nFreq)) {
        maxFreq = static_cast<std::uint64_t>(nFreq);
    }
    if (maxFreq < 1) {
        return tiles;
    }
    tiles.dofTile = modest;
    tiles.freqTile = static_cast<int>(maxFreq);
    tiles.ok = true;
    return tiles;
}

#ifdef SEMBA_DGTD_ENABLE_CUDA

// In-memory DFT. `snapshots[s]` points at 6 * nDofs doubles in component order.
// `times` are already divided by c, matching dftFromSnapshots.
// `accBudgetBytes == 0` uses free device memory. Result is divided by nSnap.
bool rcsCudaDftSnapshots(
    int nDofs,
    const double* times,
    int nSnap,
    const std::vector<double>& normFreqs,
    const double* const* snapshots,
    RcsFreqFields& ffOut,
    std::uint64_t accBudgetBytes = 0);

// Stream surface_data.bin. Applies everyNSteps and maxTime the same way as
// the host streaming DFT. Returns false when CUDA cannot hold a tile; the
// caller keeps the host DFT. Throws if the file cannot be read.
bool rcsCudaDftStreamFile(
    const std::string& surfaceBinPath,
    int nDofs,
    int everyNSteps,
    const std::optional<double>& maxTime,
    const std::vector<double>& normFreqs,
    std::vector<double>& timesOut,
    int& nSnapOut,
    RcsFreqFields& ffOut);

#endif

} // namespace maxwell
