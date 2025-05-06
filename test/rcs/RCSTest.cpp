#include <gtest/gtest.h>

#include <vector>
#include <fftw3.h>

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/RCSManager.h>
#include "components/FarField.h"

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

//TEST_F(RCSToolsTest, FunctionEval)
//{
//	Vector p({ 1.0, 1.0, 0.0 });
//	SphericalAngles angles;
//	angles.theta = M_PI_2;
//	angles.phi = 0.0;
//	Frequency freq(3e8 / physicalConstants::speedOfLight_SI);
//
//	double tol{ 1e-10 };
//	EXPECT_NEAR(0.9999905398146485235, func_exp_real_part_3D(p, freq, angles), tol);
//	EXPECT_NEAR(0.0043497449589425407, func_exp_imag_part_3D(p, freq, angles), tol);
//}

//TEST_F(RCSToolsTest, LinearFormEval)
//{
//	SphericalAngles angles;
//	angles.theta = M_PI_2;
//	angles.phi = 0.0;
//	Frequency freq(3e8 / physicalConstants::speedOfLight_SI);
//
//	Mesh mesh(1, 1, 1, Element::Type::TETRAHEDRON);
//	L2_FECollection fec(1, 3, BasisType::GaussLobatto);
//	FiniteElementSpace fes(&mesh, &fec);
//	FiniteElementSpace fes_v3(&mesh, &fec, 3);
//
//	GridFunction gf(&fes);
//	gf = 0.0;
//	gf[0] = 1.0;
//	gf[1] = 1.0;
//	gf[2] = 1.0;
//	gf[3] = 1.0;
//
//	GridFunction nodes(&fes_v3);
//	mesh.GetNodes(nodes);
//	std::vector<std::vector<double>> nodepos;
//	for (auto v = 0; v < fes.GetNDofs(); v++) {
//		nodepos.push_back({ nodes[v], nodes[v + fes.GetNDofs()], nodes[v + 2 * fes.GetNDofs()] });
//	}
//	
//	auto fc = buildFC_3D(freq, angles, true);
//	auto res{ std::make_unique<LinearForm>(&fes) };
//	res->AddBdrFaceIntegrator(new mfemExtension::FarFieldBdrFaceIntegrator(*fc.get(), X));
//	res->Assemble();
//
//}

}