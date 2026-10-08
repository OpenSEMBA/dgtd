#include "RCSSurfaceDftCuda.h"

#include "math/PhysicalConstants.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <cuda_runtime.h>

namespace maxwell {
namespace {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr int kFreqChunkLimit = 65535;
constexpr std::uint64_t kDeviceReserveBytes = 256ull << 20;

struct Dcomplex {
    double re;
    double im;
};

static_assert(sizeof(Dcomplex) == sizeof(std::complex<double>),
              "complex<double> must match the device accumulator layout");
static_assert(std::is_trivially_copyable<std::complex<double>>::value,
              "complex<double> is copied to the host as raw bytes");

void cudaCheck(cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("CUDA RCS DFT ") + what + ": " + cudaGetErrorString(err));
    }
}

struct DeviceMem {
    void* ptr{nullptr};
    explicit DeviceMem(std::size_t bytes)
    {
        if (bytes == 0) {
            return;
        }
        if (cudaMalloc(&ptr, bytes) != cudaSuccess) {
            ptr = nullptr;
        }
    }
    ~DeviceMem()
    {
        if (ptr) {
            cudaFree(ptr);
        }
    }
    DeviceMem(const DeviceMem&) = delete;
    DeviceMem& operator=(const DeviceMem&) = delete;
    bool ok() const { return ptr != nullptr; }
};

struct PinnedMem {
    void* ptr{nullptr};
    explicit PinnedMem(std::size_t bytes)
    {
        if (bytes == 0) {
            return;
        }
        if (cudaMallocHost(&ptr, bytes) != cudaSuccess) {
            ptr = nullptr;
        }
    }
    ~PinnedMem()
    {
        if (ptr) {
            cudaFreeHost(ptr);
        }
    }
    PinnedMem(const PinnedMem&) = delete;
    PinnedMem& operator=(const PinnedMem&) = delete;
    bool ok() const { return ptr != nullptr; }
};

struct Stream {
    cudaStream_t s{};
    bool live{false};
    Stream() = default;
    bool create()
    {
        if (cudaStreamCreate(&s) != cudaSuccess) {
            return false;
        }
        live = true;
        return true;
    }
    ~Stream()
    {
        if (live) {
            cudaStreamSynchronize(s);
            cudaStreamDestroy(s);
        }
    }
    Stream(const Stream&) = delete;
    Stream& operator=(const Stream&) = delete;
};

__global__ void rcsDftAccumKernel(
    Dcomplex* __restrict__ acc,
    const double* __restrict__ snap,
    const double* __restrict__ wRe,
    const double* __restrict__ wIm,
    int nDofsFull,
    int dof0,
    int nDofTile,
    int nFreqTile,
    int fBase,
    int nFreqChunk)
{
    const int dof = blockIdx.x * blockDim.x + threadIdx.x;
    const int f = blockIdx.y;
    const int c = blockIdx.z;
    if (dof >= nDofTile || f >= nFreqChunk || c >= 6) {
        return;
    }

    const double s = snap[static_cast<long long>(c) * nDofsFull + dof0 + dof];
    Dcomplex a = acc[((static_cast<long long>(c) * nFreqTile + (fBase + f)) * nDofTile) + dof];
    a.re += s * wRe[f];
    a.im += s * wIm[f];
    acc[((static_cast<long long>(c) * nFreqTile + (fBase + f)) * nDofTile) + dof] = a;
}

struct SnapSource {
    virtual ~SnapSource() = default;
    virtual void rewindRead() = 0;
    virtual bool readInto(double* dest, double& tNorm) = 0;
    virtual long long skipped() const = 0;
};

struct FileSnapSource : SnapSource {
    std::ifstream file;
    std::vector<char> iobuf;
    int nDofs{0};
    int every{1};
    std::optional<double> maxTime;
    double timeDivisor{1.0};
    std::streamoff payloadPos{0};
    long long snapIndex{0};
    long long skippedCount{0};
    bool stopped{false};

    FileSnapSource(const std::string& path, int nDofsIn, int everyIn,
                   const std::optional<double>& maxTimeIn)
        : nDofs(nDofsIn),
          every(std::max(1, everyIn)),
          maxTime(maxTimeIn),
          timeDivisor(physicalConstants::speedOfLight)
    {
        iobuf.resize(8ull << 20);
        file.rdbuf()->pubsetbuf(iobuf.data(), static_cast<std::streamsize>(iobuf.size()));
        file.open(path, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Cannot open " + path);
        }

        std::int32_t hdr[5];
        file.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
        if (!file) {
            throw std::runtime_error("Failed to read surface_data.bin header.");
        }
        if (hdr[1] != nDofs) {
            throw std::runtime_error("surface_data.bin DOF count does not match the mesh export.");
        }
        const std::uint64_t nqp = static_cast<std::uint64_t>(hdr[3]);
        const std::uint64_t sdim = static_cast<std::uint64_t>(hdr[0]);
        const std::uint64_t geoDoubles = nqp * sdim + nqp * 3ull + nqp;
        file.seekg(static_cast<std::streamoff>(geoDoubles * sizeof(double)), std::ios::cur);
        if (!file) {
            throw std::runtime_error("Failed to read surface_data.bin geometry.");
        }
        payloadPos = file.tellg();
    }

    void rewindRead() override
    {
        file.clear();
        file.seekg(payloadPos);
        snapIndex = 0;
        skippedCount = 0;
        stopped = false;
    }

    bool readInto(double* dest, double& tNorm) override
    {
        if (stopped) {
            return false;
        }
        const std::size_t payload = 6ull * static_cast<std::size_t>(nDofs) * sizeof(double);
        while (file.peek() != EOF) {
            double t = 0.0;
            file.read(reinterpret_cast<char*>(&t), sizeof(double));
            if (!file) {
                return false;
            }
            if (maxTime.has_value() && t > maxTime.value()) {
                stopped = true;
                return false;
            }
            const bool keep = (snapIndex % every) == 0;
            if (!keep) {
                file.seekg(static_cast<std::streamoff>(payload), std::ios::cur);
                if (!file) {
                    return false;
                }
                ++skippedCount;
                ++snapIndex;
                continue;
            }
            file.read(reinterpret_cast<char*>(dest), static_cast<std::streamsize>(payload));
            if (!file) {
                return false;
            }
            tNorm = t / timeDivisor;
            ++snapIndex;
            return true;
        }
        return false;
    }

    long long skipped() const override { return skippedCount; }
};

struct MemorySnapSource : SnapSource {
    const double* const* snaps{nullptr};
    const double* times{nullptr};
    int nSnap{0};
    int nDofs{0};
    int index{0};

    void rewindRead() override { index = 0; }

    bool readInto(double* dest, double& tNorm) override
    {
        if (index >= nSnap) {
            return false;
        }
        std::memcpy(dest, snaps[index], 6ull * static_cast<std::size_t>(nDofs) * sizeof(double));
        tNorm = times[index];
        ++index;
        return true;
    }

    long long skipped() const override { return 0; }
};

std::uint64_t queryAccBudgetBytes(int nDofs)
{
    std::size_t freeBytes = 0;
    std::size_t totalBytes = 0;
    if (cudaMemGetInfo(&freeBytes, &totalBytes) != cudaSuccess) {
        return 0;
    }
    const std::uint64_t snapBytes =
        6ull * static_cast<std::uint64_t>(nDofs) * sizeof(double);
    const std::uint64_t overhead = kDeviceReserveBytes + snapBytes + (1ull << 20);
    if (static_cast<std::uint64_t>(freeBytes) <= overhead) {
        return 0;
    }
    return static_cast<std::uint64_t>(freeBytes) - overhead;
}

void enqueueSnapshot(
    cudaStream_t stream,
    double* dSnap,
    const double* hSnap,
    std::size_t snapBytes,
    double* dWr,
    double* dWi,
    double* hWr,
    double* hWi,
    Dcomplex* acc,
    int nDofsFull,
    int dof0,
    int nDofTile,
    int f0,
    int nFreqTile,
    const std::vector<double>& normFreqs,
    double tNorm)
{
    cudaCheck(cudaMemcpyAsync(dSnap, hSnap, snapBytes, cudaMemcpyHostToDevice, stream),
              "snapshot copy");

    for (int fBase = 0; fBase < nFreqTile;) {
        const int nChunk = std::min(kFreqChunkLimit, nFreqTile - fBase);
        for (int i = 0; i < nChunk; ++i) {
            const double arg = 2.0 * M_PI * normFreqs[f0 + fBase + i] * tNorm;
            hWr[i] = std::cos(arg);
            hWi[i] = -std::sin(arg);
        }
        const std::size_t wBytes = static_cast<std::size_t>(nChunk) * sizeof(double);
        cudaCheck(cudaMemcpyAsync(dWr, hWr, wBytes, cudaMemcpyHostToDevice, stream),
                  "weight copy");
        cudaCheck(cudaMemcpyAsync(dWi, hWi, wBytes, cudaMemcpyHostToDevice, stream),
                  "weight copy");

        const dim3 block(256);
        const dim3 grid((nDofTile + 255) / 256, nChunk, 6);
        rcsDftAccumKernel<<<grid, block, 0, stream>>>(
            acc, dSnap, dWr, dWi,
            nDofsFull, dof0, nDofTile, nFreqTile, fBase, nChunk);
        cudaCheck(cudaGetLastError(), "DFT kernel launch");

        fBase += nChunk;
        if (fBase < nFreqTile) {
            cudaCheck(cudaStreamSynchronize(stream), "frequency-chunk sync");
        }
    }
}

bool runTiledDft(
    SnapSource& source,
    int nDofs,
    int everyNSteps,
    const std::vector<double>& normFreqs,
    std::uint64_t accBudgetBytes,
    bool printPlan,
    std::vector<double>& timesOut,
    int& nSnapOut,
    RcsFreqFields& ffOut)
{
    if (nDofs <= 0 || normFreqs.empty() ||
        normFreqs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return false;
    }
    const int nFreq = static_cast<int>(normFreqs.size());
    const std::uint64_t budget =
        accBudgetBytes == 0 ? queryAccBudgetBytes(nDofs) : accBudgetBytes;
    const RcsCudaDftTiles tiles = chooseRcsCudaDftTiles(nFreq, nDofs, budget);
    if (!tiles.ok) {
        if (printPlan) {
            std::cout << "    CUDA DFT   : accumulator does not fit in free VRAM; "
                      << "using host DFT.\n";
        }
        return false;
    }

    const std::size_t snapBytes =
        6ull * static_cast<std::size_t>(nDofs) * sizeof(double);
    const int freqChunkCap = std::min(tiles.freqTile, kFreqChunkLimit);
    const std::size_t wBytes = static_cast<std::size_t>(freqChunkCap) * sizeof(double);
    const std::size_t accBytes =
        6ull * static_cast<std::size_t>(tiles.freqTile) *
        static_cast<std::size_t>(tiles.dofTile) * sizeof(Dcomplex);

    PinnedMem pin0(snapBytes);
    PinnedMem pin1(snapBytes);
    PinnedMem hWr(wBytes);
    PinnedMem hWi(wBytes);
    DeviceMem dSnap(snapBytes);
    DeviceMem dWr(wBytes);
    DeviceMem dWi(wBytes);
    DeviceMem acc(accBytes);
    Stream stream;
    if (!pin0.ok() || !pin1.ok() || !hWr.ok() || !hWi.ok() ||
        !dSnap.ok() || !dWr.ok() || !dWi.ok() || !acc.ok() || !stream.create()) {
        if (printPlan) {
            std::cout << "    CUDA DFT   : device allocation failed; using host DFT.\n";
        }
        return false;
    }

    const int passes =
        rcsCudaDftPassCount(nDofs, tiles.dofTile) *
        rcsCudaDftPassCount(nFreq, tiles.freqTile);
    if (printPlan) {
        std::cout << "    CUDA DFT   : tile " << tiles.dofTile << " DOFs x "
                  << tiles.freqTile << " frequencies, " << passes
                  << " file pass(es)\n" << std::flush;
    }

    RcsFreqFields ff(6, std::vector<std::vector<std::complex<double>>>(
        static_cast<std::size_t>(nFreq),
        std::vector<std::complex<double>>(static_cast<std::size_t>(nDofs), {0.0, 0.0})));

    double* pins[2] = {
        static_cast<double*>(pin0.ptr),
        static_cast<double*>(pin1.ptr)
    };
    auto* hWrPtr = static_cast<double*>(hWr.ptr);
    auto* hWiPtr = static_cast<double*>(hWi.ptr);
    auto* dSnapPtr = static_cast<double*>(dSnap.ptr);
    auto* dWrPtr = static_cast<double*>(dWr.ptr);
    auto* dWiPtr = static_cast<double*>(dWi.ptr);
    auto* accPtr = static_cast<Dcomplex*>(acc.ptr);

    int nSnap = -1;
    std::vector<double> times;
    const int stride = std::max(1, everyNSteps);

    for (int dof0 = 0; dof0 < nDofs; dof0 += tiles.dofTile) {
        const int nDofTile = std::min(tiles.dofTile, nDofs - dof0);
        for (int freq0 = 0; freq0 < nFreq; freq0 += tiles.freqTile) {
            const int nFreqTile = std::min(tiles.freqTile, nFreq - freq0);
            const std::size_t usedBytes =
                6ull * static_cast<std::size_t>(nFreqTile) *
                static_cast<std::size_t>(nDofTile) * sizeof(Dcomplex);
            cudaCheck(cudaMemsetAsync(accPtr, 0, usedBytes, stream.s), "accumulator clear");

            source.rewindRead();
            int count = 0;
            std::vector<double> passTimes;
            double tNorm = 0.0;
            int slot = 0;
            bool more = source.readInto(pins[0], tNorm);
            while (more) {
                enqueueSnapshot(
                    stream.s, dSnapPtr, pins[slot], snapBytes,
                    dWrPtr, dWiPtr, hWrPtr, hWiPtr, accPtr,
                    nDofs, dof0, nDofTile, freq0, nFreqTile, normFreqs, tNorm);
                if (nSnap < 0) {
                    passTimes.push_back(tNorm);
                }
                ++count;
                const int nextSlot = 1 - slot;
                double tNext = 0.0;
                more = source.readInto(pins[nextSlot], tNext);
                cudaCheck(cudaStreamSynchronize(stream.s), "snapshot sync");
                if (!more) {
                    break;
                }
                slot = nextSlot;
                tNorm = tNext;
            }

            if (count == 0) {
                cudaCheck(cudaStreamSynchronize(stream.s), "empty-stream sync");
                throw std::runtime_error("No snapshots found while streaming surface_data.bin.");
            }
            if (nSnap < 0) {
                nSnap = count;
                times = std::move(passTimes);
                if (printPlan && stride > 1) {
                    std::cout << "    every_n_steps=" << stride
                              << ": kept " << nSnap << ", skipped " << source.skipped()
                              << " snapshots while streaming\n";
                }
            } else if (count != nSnap) {
                throw std::runtime_error("CUDA RCS DFT passes saw different snapshot counts.");
            }

            cudaCheck(cudaStreamSynchronize(stream.s), "pass sync");
            std::vector<std::complex<double>> host(
                6ull * static_cast<std::size_t>(nFreqTile) * static_cast<std::size_t>(nDofTile));
            cudaCheck(cudaMemcpy(host.data(), accPtr,
                                 host.size() * sizeof(std::complex<double>),
                                 cudaMemcpyDeviceToHost),
                      "accumulator copy");
            for (int c = 0; c < 6; ++c) {
                for (int lf = 0; lf < nFreqTile; ++lf) {
                    auto& row = ff[c][freq0 + lf];
                    const std::complex<double>* src =
                        host.data() + (static_cast<std::size_t>(c) * nFreqTile + lf) * nDofTile;
                    for (int ld = 0; ld < nDofTile; ++ld) {
                        row[dof0 + ld] = src[ld];
                    }
                }
            }
        }
    }

    const double invN = 1.0 / static_cast<double>(nSnap);
    for (auto& comp : ff) {
        for (auto& row : comp) {
            for (auto& v : row) {
                v *= invN;
            }
        }
    }

    timesOut = std::move(times);
    nSnapOut = nSnap;
    ffOut = std::move(ff);
    return true;
}

} // namespace

bool rcsCudaDftSnapshots(
    int nDofs,
    const double* times,
    int nSnap,
    const std::vector<double>& normFreqs,
    const double* const* snapshots,
    RcsFreqFields& ffOut,
    std::uint64_t accBudgetBytes)
{
    if (nSnap <= 0 || times == nullptr || snapshots == nullptr) {
        return false;
    }
    MemorySnapSource source;
    source.snaps = snapshots;
    source.times = times;
    source.nSnap = nSnap;
    source.nDofs = nDofs;
    std::vector<double> timesOut;
    int nSnapOut = 0;
    return runTiledDft(source, nDofs, 1, normFreqs, accBudgetBytes, false,
                       timesOut, nSnapOut, ffOut);
}

bool rcsCudaDftStreamFile(
    const std::string& surfaceBinPath,
    int nDofs,
    int everyNSteps,
    const std::optional<double>& maxTime,
    const std::vector<double>& normFreqs,
    std::vector<double>& timesOut,
    int& nSnapOut,
    RcsFreqFields& ffOut)
{
    FileSnapSource source(surfaceBinPath, nDofs, everyNSteps, maxTime);
    return runTiledDft(source, nDofs, everyNSteps, normFreqs, 0, true,
                       timesOut, nSnapOut, ffOut);
}

} // namespace maxwell
