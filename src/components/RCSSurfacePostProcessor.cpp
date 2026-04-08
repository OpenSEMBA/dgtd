#include "RCSSurfacePostProcessor.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <regex>
#include <stdexcept>

#include <omp.h>

#include "driver/driver.h"
#include "math/PhysicalConstants.h"
#include "mfemExtension/LinearIntegrators.h"

namespace maxwell {

using namespace mfem;

// ------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------

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

// Spherical coordinate unit vectors.
static std::array<double, 3> thetaHat(const SphericalAngles& a)
{
    return { std::cos(a.theta)*std::cos(a.phi),
             std::cos(a.theta)*std::sin(a.phi),
            -std::sin(a.theta) };
}

static std::array<double, 3> phiHat(const SphericalAngles& a)
{
    return { -std::sin(a.phi), std::cos(a.phi), 0.0 };
}

static std::array<double, 3> rHat(const SphericalAngles& a)
{
    return { std::sin(a.theta)*std::cos(a.phi),
             std::sin(a.theta)*std::sin(a.phi),
             std::cos(a.theta) };
}

template <typename T>
static T dot3(const std::array<double,3>& a, const std::array<T,3>& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
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

    int32_t hdr[4];
    f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    rd.geometry.spaceDimension = hdr[0];
    rd.geometry.numDofs        = hdr[1];
    rd.geometry.numBdrElements = hdr[2];
    rd.geometry.numQuadPoints  = hdr[3];

    const int nqp = rd.geometry.numQuadPoints;
    const int sdim = rd.geometry.spaceDimension;

    rd.geometry.positions.resize(nqp * sdim);
    rd.geometry.normals.resize(nqp * 3);
    rd.geometry.weights.resize(nqp);

    f.read(reinterpret_cast<char*>(rd.geometry.positions.data()),
           nqp * sdim * sizeof(double));
    f.read(reinterpret_cast<char*>(rd.geometry.normals.data()),
           nqp * 3 * sizeof(double));
    f.read(reinterpret_cast<char*>(rd.geometry.weights.data()),
           nqp * sizeof(double));

    const int ndofs = rd.geometry.numDofs;
    while (f.peek() != EOF) {
        SurfaceSnapshot snap;
        f.read(reinterpret_cast<char*>(&snap.time), sizeof(double));
        if (!f) break;

        for (auto* vec : {&snap.Ex, &snap.Ey, &snap.Ez,
                          &snap.Hx, &snap.Hy, &snap.Hz}) {
            vec->resize(ndofs);
            f.read(reinterpret_cast<char*>(vec->data()), ndofs * sizeof(double));
        }
        rd.snapshots.push_back(std::move(snap));
    }
    return rd;
}

// ------------------------------------------------------------------
// Plane-wave helpers
// ------------------------------------------------------------------

PlaneWaveData RCSSurfacePostProcessor::extractPlaneWaveData(
    const std::string& jsonPath) const
{
    return buildPlaneWaveData(driver::parseJSONfile(jsonPath));
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

    // Rescale frequencies from Hz to normalised (f / c_SI).
    std::vector<double> normFreqs(frequencies.size());
    for (size_t i = 0; i < frequencies.size(); ++i)
        normFreqs[i] = frequencies[i] / physicalConstants::speedOfLight_SI;

    const int nFreq = static_cast<int>(normFreqs.size());

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

    for (const auto& rp : rankPaths) {
        auto rd = readRankData(rp);
        spaceDim = rd.geometry.spaceDimension;
        const int nDofs = rd.geometry.numDofs;

        // Filter snapshots based on maxTime if provided.
        std::vector<SurfaceSnapshot> snapshots = rd.snapshots;
        if (maxTime_.has_value()) {
            snapshots.erase(
                std::remove_if(snapshots.begin(), snapshots.end(),
                    [this](const SurfaceSnapshot& s) { return s.time > maxTime_.value(); }),
                snapshots.end());
            if (snapshots.empty()) {
                throw std::runtime_error("No snapshots found within the specified maxTime.");
            }
        }
        const int nSnap = static_cast<int>(snapshots.size());

        if (spaceDim == 2)
            for (const auto& a : angles)
                if (std::abs(a.theta - M_PI_2) > 1e-8)
                    throw std::runtime_error("2D RCS requires theta = pi/2.");

        // Time vector (normalised).
        std::vector<double> times(nSnap);
        for (int i = 0; i < nSnap; ++i)
            times[i] = snapshots[i].time / physicalConstants::speedOfLight;

        if (firstRank) {
            pw = extractPlaneWaveData(jsonPath);
            incidentPower = computeIncidentPowerSpectrum(pw, times, normFreqs);
            firstRank = false;
        }

        // --- DFT all 6 field components ---
        // freqFields[comp][freq] is a ComplexVector of size nDofs.
        // comp: 0=Ex 1=Ey 2=Ez 3=Hx 4=Hy 5=Hz
        using CVec = std::vector<std::complex<double>>;
        std::vector<std::vector<CVec>> ff(6,
            std::vector<CVec>(nFreq, CVec(nDofs, {0,0})));

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

        // --- Load mesh and build FES for this rank ---
        auto mesh = Mesh::LoadFromFile(rp + "/mesh", 1, 0);
        auto pmesh = ParMesh(MPI_COMM_WORLD, mesh);
        int order = determineFECOrder(pmesh, nDofs);
        DG_FECollection fec(order, pmesh.Dimension());
        ParFiniteElementSpace fes(&pmesh, &fec);

        // --- For each frequency and angle, compute far-field potentials ---
        for (int fi = 0; fi < nFreq; ++fi) {
            const double freq = normFreqs[fi];
            const double k = 2.0 * M_PI * freq;

            for (const auto& ang : angles) {
                // Build phase-term function coefficients.
                std::unique_ptr<FunctionCoefficient> fcR, fcI;
                if (spaceDim == 2) {
                    fcR = buildFC_2D(freq, ang, true);
                    fcI = buildFC_2D(freq, ang, false);
                } else {
                    fcR = buildFC_3D(freq, ang, true);
                    fcI = buildFC_3D(freq, ang, false);
                }

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
                auto th = thetaHat(ang);
                auto ph = phiHat(ang);

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
    }

    // After all ranks processed, compute far-field power from coherently-summed amplitudes.
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

    // Compute RCS from accumulated far-field data.
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

    // Write merged output files.
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
}

// ------------------------------------------------------------------
// Constructor
// ------------------------------------------------------------------

RCSSurfacePostProcessor::RCSSurfacePostProcessor(
    const std::string& dataPath,
    const std::string& jsonPath,
    std::vector<Frequency>& frequencies,
    const std::vector<SphericalAngles>& angles,
    const std::optional<double>& maxTime)
    : maxTime_(maxTime)
{
    computeAndWriteResults(dataPath, jsonPath, frequencies, angles);
}

} // namespace maxwell
