#include <gtest/gtest.h>

#include <vector>
#include <fftw3.h>

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/RCSExporter.h>
#include <components/RCSManager.h>
#include <components/RCSManager.cpp>

namespace maxwell {

using namespace mfem;

class RCSTest : public ::testing::Test {
};

TEST_F(RCSTest, circleTest)
{
	const double dt = 5e-3;
	const int steps = 100000;

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {M_PI * 3.0 / 4.0, 0.0} };

	RCSManager rcs("NearToFarFieldExports/circle", "2D_NTFF_Circle", dt, steps, angles);
}

TEST_F(RCSTest, sphereTest)
{
	const double dt = 5e-3;
	const int steps = 100000;

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/sphere", "3D_RCS", dt, steps, angles);
}

TEST_F(RCSTest, DiscreteFourierTransform)
{
	double in[3] = { 1.0, 2.0, 3.0 };
	Vector field(in);
	std::vector<double> times({ 5e-3, 10e-3 });
	const int N = 3;
	fftw_complex* out;
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 3);

	fftw_plan p;

	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

	fftw_execute(p);

	std::vector<double> frequencies;
	std::vector<double> fftw_mag;
	for (int i = 0; i < N / 2 + 1; ++i) {
		frequencies.push_back(i * (1.0 / (times[1] - times[0])) / N);
		fftw_mag.push_back(sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]));
	}

	auto dft_0{ calculateDFT(field, frequencies, times[0])};
	auto dft_1{ calculateDFT(field, frequencies, times[1])};
	std::vector<double> dft_mag;
	dft_mag.push_back(sqrt(dft_0[0][0].real() * dft_0[0][0].real() + dft_0[0][0].imag() * dft_0[0][0].imag() + dft_0[1][0].real() * dft_0[1][0].real() + dft_0[1][0].imag() * dft_0[1][0].imag()));


	fftw_destroy_plan(p);
	fftw_free(out);
}

}