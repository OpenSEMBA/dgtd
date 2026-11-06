#include <gtest/gtest.h>

#include <vector>
#include <fftw3.h>

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <cmath>
#include <complex>
#include <cstdint>

#include <components/RCSManager.h>
#include "components/FarField.h"
#include "components/RCSSurfaceBinaryWriter.h"
#include "components/RCSSurfaceDftCuda.h"
#include "components/RCSSurfacePostProcessor.h"
#include "driver/driver.h"

#include "TestUtils.h"

namespace maxwell {

using namespace mfem;

class RCSToolsTest : public ::testing::Test{
};

TEST_F(RCSToolsTest, DiscreteFourierTransform)
{
	const int N = 3;
	double in[N] = { 1.0, 2.0, 3.0 };
	Vector field(in);
	std::vector<double> times({ 5e-3, 10e-3 });
	fftw_complex* out;
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	fftw_plan p;

	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

	fftw_execute(p);

	std::vector<double> frequencies;
	std::vector<double> fftw_mag;
	for (int i = 0; i < N / 2 + 1; ++i) {
		frequencies.push_back(i * (1.0 / (times[1] - times[0])) / N);
		fftw_mag.push_back(sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]));
	}

	auto dft_0{ calculateDFT(field, frequencies, times[0]) };
	auto dft_1{ calculateDFT(field, frequencies, times[1]) };
	std::vector<double> dft_mag;
	dft_mag.push_back(sqrt(dft_0[0][0].real() * dft_0[0][0].real() + dft_0[0][0].imag() * dft_0[0][0].imag() + dft_0[1][0].real() * dft_0[1][0].real() + dft_0[1][0].imag() * dft_0[1][0].imag()));


	fftw_destroy_plan(p);
	fftw_free(out);
}

TEST_F(RCSToolsTest, PsiAngle) 
{
	SphericalAngles angleDirX;
	angleDirX.theta = M_PI_2;
	angleDirX.phi = 0.0;
	Vector vecDirX({ 1.0, 0.0, 0.0 });
	Vector vecDirY({ 0.0, 1.0, 0.0 });
	Vector vecDirZ({ 0.0, 0.0, 1.0 });
	Vector vecDirXY({ 1.0, 1.0, 0.0 });
	EXPECT_EQ(0.0, calcPsiAngle3D(vecDirX, angleDirX));
	double tol{ 1e-8 };
	EXPECT_NEAR(M_PI_2, calcPsiAngle3D(vecDirY, angleDirX), tol);
	EXPECT_NEAR(M_PI_2, calcPsiAngle3D(vecDirZ, angleDirX), tol);
	EXPECT_NEAR(M_PI_4, calcPsiAngle3D(vecDirXY, angleDirX), tol);
}

TEST_F(RCSToolsTest, FunctionEval)
{
	Vector p({ 1.0, 1.0, 0.0 });
	SphericalAngles angles;
	angles.theta = M_PI_2;
	angles.phi = 0.0;
	Frequency freq(3e8 / physicalConstants::speedOfLight_SI);

	double tol{ 1e-10 };
	EXPECT_NEAR(0.9999905398146485235, evalFuncExpPart(p, freq, angles, 3, true), tol);
	EXPECT_NEAR(0.0043497449589425407, evalFuncExpPart(p, freq, angles, 3, false), tol);
}

TEST_F(RCSToolsTest, LinearFormEval)
{
	SphericalAngles angles;
	angles.theta = M_PI_2;
	angles.phi = 0.0;
	Frequency freq(3e8 / physicalConstants::speedOfLight_SI);

	Mesh mesh = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::TETRAHEDRON);
	L2_FECollection fec(1, 3, BasisType::GaussLobatto);
	FiniteElementSpace fes(&mesh, &fec);
	FiniteElementSpace fes_v3(&mesh, &fec, 3);

	GridFunction gf(&fes);
	gf = 0.0;
	gf[0] = 1.0;
	gf[1] = 1.0;
	gf[2] = 1.0;
	gf[3] = 1.0;

	GridFunction nodes(&fes_v3);
	mesh.GetNodes(nodes);
	std::vector<std::vector<double>> nodepos;
	for (auto v = 0; v < fes.GetNDofs(); v++) {
		nodepos.push_back({ nodes[v], nodes[v + fes.GetNDofs()], nodes[v + 2 * fes.GetNDofs()] });
	}
	
	auto fc = buildFC(3, freq, angles, true);
	auto res{ std::make_unique<LinearForm>(&fes) };
	res->AddBdrFaceIntegrator(new mfemExtension::FarFieldBdrFaceIntegrator(*fc.get(), X));
	res->Assemble();

}

TEST_F(RCSToolsTest, CudaDftTileChoice)
{
	const auto bytes = rcsCudaDftAccBytesPerFreqDof();
	EXPECT_EQ(96ull, bytes);

	const auto full = chooseRcsCudaDftTiles(4, 1000, bytes * 4ull * 1000ull);
	EXPECT_TRUE(full.ok);
	EXPECT_EQ(1000, full.dofTile);
	EXPECT_EQ(4, full.freqTile);
	EXPECT_EQ(1, rcsCudaDftPassCount(1000, full.dofTile));

	const auto dofTiled = chooseRcsCudaDftTiles(4, 1000, bytes * 4ull * 500ull);
	EXPECT_TRUE(dofTiled.ok);
	EXPECT_EQ(4, dofTiled.freqTile);
	EXPECT_EQ(384, dofTiled.dofTile);
	EXPECT_EQ(3, rcsCudaDftPassCount(1000, dofTiled.dofTile));

	const auto freqTiled = chooseRcsCudaDftTiles(100000, 100000, bytes * 10ull * 4096ull);
	EXPECT_TRUE(freqTiled.ok);
	EXPECT_EQ(4096, freqTiled.dofTile);
	EXPECT_EQ(10, freqTiled.freqTile);

	const auto tooSmall = chooseRcsCudaDftTiles(3, 5, bytes - 1);
	EXPECT_FALSE(tooSmall.ok);
}

#ifdef SEMBA_DGTD_ENABLE_CUDA
namespace {

RcsFreqFields referenceDft(
	int nSnap, int nFreq, int nDofs,
	const std::vector<double>& times,
	const std::vector<double>& freqs,
	const std::vector<std::vector<double>>& snaps)
{
	RcsFreqFields ff(6, std::vector<std::vector<std::complex<double>>>(
		static_cast<std::size_t>(nFreq),
		std::vector<std::complex<double>>(static_cast<std::size_t>(nDofs), {0.0, 0.0})));
	for (int s = 0; s < nSnap; ++s) {
		for (int fi = 0; fi < nFreq; ++fi) {
			const double arg = 2.0 * M_PI * freqs[fi] * times[s];
			const std::complex<double> w(std::cos(arg), -std::sin(arg));
			for (int c = 0; c < 6; ++c) {
				for (int v = 0; v < nDofs; ++v) {
					ff[c][fi][v] += snaps[s][c * nDofs + v] * w;
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
	return ff;
}

void expectDftMatch(const RcsFreqFields& got, const RcsFreqFields& ref)
{
	ASSERT_EQ(ref.size(), got.size());
	for (std::size_t c = 0; c < ref.size(); ++c) {
		ASSERT_EQ(ref[c].size(), got[c].size());
		for (std::size_t fi = 0; fi < ref[c].size(); ++fi) {
			ASSERT_EQ(ref[c][fi].size(), got[c][fi].size());
			for (std::size_t v = 0; v < ref[c][fi].size(); ++v) {
				EXPECT_NEAR(ref[c][fi][v].real(), got[c][fi][v].real(), 1e-12);
				EXPECT_NEAR(ref[c][fi][v].imag(), got[c][fi][v].imag(), 1e-12);
			}
		}
	}
}

} // namespace

TEST_F(RCSToolsTest, CudaDiscreteFourierTransform)
{
	if (!Device::Allows(Backend::CUDA)) {
		GTEST_SKIP() << "CUDA backend not active.";
	}

	const int nSnap = 5;
	const int nFreq = 4;
	const int nDofs = 7;
	const std::vector<double> times{0.0, 0.05, 0.10, 0.15, 0.20};
	const std::vector<double> freqs{0.25, 1.0, 2.5, 4.0};
	std::vector<std::vector<double>> storage(
		nSnap, std::vector<double>(6 * nDofs));
	std::vector<const double*> ptrs(nSnap);
	for (int s = 0; s < nSnap; ++s) {
		ptrs[s] = storage[s].data();
		for (int c = 0; c < 6; ++c) {
			for (int v = 0; v < nDofs; ++v) {
				storage[s][c * nDofs + v] =
					0.1 * (s + 1) + 0.01 * (c + 1) + 0.001 * v;
			}
		}
	}
	const auto ref = referenceDft(nSnap, nFreq, nDofs, times, freqs, storage);

	RcsFreqFields full;
	ASSERT_TRUE(rcsCudaDftSnapshots(
		nDofs, times.data(), nSnap, freqs, ptrs.data(), full, 0));
	expectDftMatch(full, ref);

	const auto bytes = rcsCudaDftAccBytesPerFreqDof();
	RcsFreqFields tiled;
	ASSERT_TRUE(rcsCudaDftSnapshots(
		nDofs, times.data(), nSnap, freqs, ptrs.data(), tiled,
		bytes * static_cast<std::uint64_t>(nFreq) * 3ull));
	expectDftMatch(tiled, ref);
}
#endif

TEST_F(RCSToolsTest, RCSSurfacePostProcessor_LoadAllMatchesRamGatedStreamingAndFilters)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Synthetic single-rank RCSSurface dump check runs on 1 rank.";
	}

	const std::string case_path = maxwellCase("2D_RCS_Circle_G1");
	auto case_data = driver::parseJSONfile(case_path);
	auto model = driver::buildModel(case_data, case_path, true);

	DG_FECollection fec(1, 2, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&model.getMesh(), &fec);
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);

	const std::string tmp_dir = "rcs_stream_test_tmp";
	const std::string rank_dir = tmp_dir + "/rank0";
	std::filesystem::remove_all(tmp_dir);

	{
		RCSSurfaceBinaryWriter writer({2, 3, 4, 5}, &fec, fes, fields, rank_dir);
		const int ndofs = fes.GetNDofs();
		for (int step = 0; step < 8; ++step) {
			const double t = 0.05 * static_cast<double>(step);
			for (int i = 0; i < ndofs; ++i) {
				fields.get(E, Y)[i] = std::sin(2.0 * M_PI * (t + 0.01 * i));
				fields.get(H, Z)[i] = std::cos(2.0 * M_PI * (t + 0.01 * i));
			}
			writer.writeSnapshot(t);
		}
		writer.flush();
	}

	std::vector<Frequency> freqs = {3.0e8, 6.0e8};
	SphericalAngles ang{M_PI_2, M_PI};
	std::vector<SphericalAngles> angles = {ang};

	// 1. Load-all (ram_gate = 10.0 GiB) vs RAM-gated streaming DFT (ram_gate = 1e-9 GiB)
	RCSSurfacePostProcessor pp_load_all(
		tmp_dir, case_path, freqs, angles, std::nullopt, 1, 10.0);
	RCSSurfacePostProcessor pp_stream(
		tmp_dir, case_path, freqs, angles, std::nullopt, 1, 1.0e-9);

	for (double f_hz : freqs) {
		const double f_norm = f_hz / physicalConstants::speedOfLight_SI;
		EXPECT_GT(pp_load_all.getRCS(ang, f_norm), 0.0);
		EXPECT_NEAR(
			pp_load_all.getRCS(ang, f_norm),
			pp_stream.getRCS(ang, f_norm),
			1e-12);
		EXPECT_NEAR(
			pp_load_all.getFarField(ang, f_norm),
			pp_stream.getFarField(ang, f_norm),
			1e-12);
	}

	// 2. Load-all vs streaming with max_time and every_n_steps filters
	RCSSurfacePostProcessor pp_filtered_load(
		tmp_dir, case_path, freqs, angles, 0.22, 2, 10.0);
	RCSSurfacePostProcessor pp_filtered_stream(
		tmp_dir, case_path, freqs, angles, 0.22, 2, 1.0e-9);

	for (double f_hz : freqs) {
		const double f_norm = f_hz / physicalConstants::speedOfLight_SI;
		EXPECT_NEAR(
			pp_filtered_load.getRCS(ang, f_norm),
			pp_filtered_stream.getRCS(ang, f_norm),
			1e-12);
	}

	// 3. 2D RCS must reject theta != pi/2
	std::vector<SphericalAngles> bad_angles = {{0.5, 0.0}};
	EXPECT_THROW(
		RCSSurfacePostProcessor(tmp_dir, case_path, freqs, bad_angles, std::nullopt, 1, 10.0),
		std::runtime_error);

	std::filesystem::remove_all(tmp_dir);
}

}