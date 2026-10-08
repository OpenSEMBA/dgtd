#include "math/Function.h"
#include <gtest/gtest.h>
#include <cmath>

class FunctionTest : public ::testing::Test {
};

using namespace maxwell;
using namespace mfem;

TEST_F(FunctionTest, DerivGaussianTest) {
	
	auto dgauss{ DerivGaussDipole(0.1, 0.2, 0.8) };
	auto dt{ 1e-3 };
	auto final_t{ 3.0 };
	auto steps = int(final_t / dt);
	std::vector<double> time_vector(steps);
	for (auto t{ 0 }; t < time_vector.size(); t++) {
		time_vector[t] = t * dt;
	}

	Position p({ 0.0,0.0,0.314158 });
	std::vector<double> resx(steps), resy(steps), resz(steps);
	for (auto t{ 0.0 }; t < time_vector.size(); t++) {
		resx[t] = dgauss.eval(p, time_vector[t], FieldType::E, X);
		resy[t] = dgauss.eval(p, time_vector[t], FieldType::E, Y);
		resz[t] = dgauss.eval(p, time_vector[t], FieldType::E, Z);
	}
}

TEST_F(FunctionTest, Gaussian_1D_peakAtMean)
{
	Position mean({0.0});
	Gaussian g(1.0, mean, 1);

	Position atMean({0.0});
	EXPECT_DOUBLE_EQ(1.0, g.eval(atMean));
	EXPECT_EQ(1, g.dimension());
}

TEST_F(FunctionTest, Gaussian_1D_knownValue)
{
	Position mean({0.0});
	Gaussian g(1.0, mean, 1);

	// At x = spread, value = exp(-0.5)
	Position atSpread({1.0});
	EXPECT_NEAR(std::exp(-0.5), g.eval(atSpread), 1e-15);
}

TEST_F(FunctionTest, Gaussian_1D_symmetry)
{
	Position mean({2.0});
	Gaussian g(0.5, mean, 1);

	Position left({1.5});
	Position right({2.5});
	EXPECT_DOUBLE_EQ(g.eval(left), g.eval(right));
}

TEST_F(FunctionTest, Gaussian_2D_peakAtMean)
{
	Position mean({0.0, 0.0});
	Gaussian g(1.0, mean, 2);

	EXPECT_EQ(2, g.dimension());
	Position atOrigin({0.0, 0.0});
	EXPECT_DOUBLE_EQ(1.0, g.eval(atOrigin));
}

TEST_F(FunctionTest, Gaussian_3D_peakAtMean)
{
	Position mean({0.0, 0.0, 0.0});
	Gaussian g(1.0, mean, 3);

	EXPECT_EQ(3, g.dimension());
	Position atOrigin({0.0, 0.0, 0.0});
	EXPECT_DOUBLE_EQ(1.0, g.eval(atOrigin));
}

TEST_F(FunctionTest, ModulatedGaussian_atMean)
{
	Position mean({0.5});
	ModulatedGaussian mg(0.2, mean, 3.0, 1);

	// At x = mean: exp(0) * cos(0) = 1.0
	Position atMean({0.5});
	EXPECT_DOUBLE_EQ(1.0, mg.eval(atMean));
	EXPECT_EQ(1, mg.dimension());
}

TEST_F(FunctionTest, ModulatedGaussian_zeroCrossing)
{
	Position mean({0.0});
	double spread = 10.0; // wide envelope
	double freq = 1.0;
	ModulatedGaussian mg(spread, mean, freq, 1);

	// At x = 1/(4*freq) = 0.25, cos(2*pi*freq*x) = cos(pi/2) = 0
	Position quarterWavelength({0.25});
	EXPECT_NEAR(0.0, mg.eval(quarterWavelength), 1e-10);
}

TEST_F(FunctionTest, SinusoidalMode_boundaryZeros_1D)
{
	SinusoidalMode mode({1});

	Position atZero({0.0});
	Position atOne({1.0});
	// sin(pi*0) = 0, sin(pi*1) ~ 0
	EXPECT_NEAR(0.0, mode.eval(atZero), 1e-15);
	EXPECT_NEAR(0.0, mode.eval(atOne),  1e-15);
}

TEST_F(FunctionTest, SinusoidalMode_peakAtHalf)
{
	SinusoidalMode mode({1});

	Position atHalf({0.5});
	// sin(pi*0.5) = 1.0
	EXPECT_NEAR(1.0, mode.eval(atHalf), 1e-15);
}

TEST_F(FunctionTest, SinusoidalMode_2D)
{
	SinusoidalMode mode({1, 2});
	EXPECT_EQ(2, mode.dimension());

	Position p({0.5, 0.25});
	// sin(pi*0.5) * sin(2*pi*0.25) = 1.0 * 1.0 = 1.0
	EXPECT_NEAR(1.0, mode.eval(p), 1e-15);
}

TEST_F(FunctionTest, SinusoidalMode_higherModes)
{
	SinusoidalMode mode({2});

	// sin(2*pi*0.25) = sin(pi/2) = 1.0
	Position p({0.25});
	EXPECT_NEAR(1.0, mode.eval(p), 1e-15);

	// sin(2*pi*0.5) = sin(pi) ~ 0
	Position half({0.5});
	EXPECT_NEAR(0.0, mode.eval(half), 1e-15);
}

TEST_F(FunctionTest, GaussianSpread_power1eMatchesHz)
{
	const double f_1e = 3.0e8;
	const double sigma = gaussianSpreadForPower1e(f_1e);
	const double omega = 2.0 * M_PI * f_1e / physicalConstants::speedOfLight_SI;
	EXPECT_NEAR(1.0 / std::exp(1.0), std::exp(-(omega * sigma) * (omega * sigma)), 1e-12);
	EXPECT_NEAR(0.15, gaussianSpreadForPower1e(
		physicalConstants::speedOfLight_SI / (2.0 * M_PI * 0.15)), 1e-12);
}

TEST_F(FunctionTest, GaussianSpread_minus20dB)
{
	const double f_max = 0.01;
	const double sigma = gaussianSpreadForDbCut(f_max, -20.0);
	const double pi = std::acos(-1.0);
	const double ratio = std::exp(-2.0 * pi * pi * sigma * sigma * f_max * f_max);
	EXPECT_NEAR(0.1, ratio, 1e-12);
}

TEST_F(FunctionTest, GaussianSpread_minus30dB)
{
	const double f_max = 2.0;
	const double sigma = gaussianSpreadForDbCut(f_max, -30.0);
	const double pi = std::acos(-1.0);
	const double ratio = std::exp(-2.0 * pi * pi * sigma * sigma * f_max * f_max);
	EXPECT_NEAR(std::pow(10.0, -1.5), ratio, 1e-12);
}

TEST_F(FunctionTest, GaussianTimeSignal_derivativeAccumulatesThePulse)
{
	const double magnitude = 1.0;
	const double spread = 0.5;
	const double t0 = 2.0;
	EXPECT_NEAR(magnitude, gaussianTimeSignal(t0, t0, spread, magnitude, false), 1e-12);
	EXPECT_NEAR(0.0, gaussianTimeSignal(t0, t0, spread, magnitude, true), 1e-12);
	const double flank = t0 - spread;
	EXPECT_NEAR(
		magnitude * std::exp(-0.5) / spread,
		gaussianTimeSignal(flank, t0, spread, magnitude, true),
		1e-12);

	double integral = 0.0;
	const double dt = spread / 200.0;
	for (double t = t0 - 12.0 * spread; t < t0 + dt * 0.5; t += dt) {
		integral += gaussianTimeSignal(t, t0, spread, magnitude, true) * dt;
	}
	EXPECT_NEAR(magnitude, integral, 1e-3);
}

TEST_F(FunctionTest, CoaxialMode_radialElectricAndAzimuthalMagnetic)
{
	Vector center(3);
	center = 0.0;
	Vector axis({1.0, 0.0, 0.0});
	const double a = 0.05;
	const double b = 0.10;
	const double v0 = 1.0;
	CoaxialMode mode(v0, 1.0, 0.0, false, center, axis, a, b, 113, 112);
	Vector at({0.0, 0.075, 0.0});
	const double e_rho = v0 / (0.075 * std::log(b / a));
	EXPECT_NEAR(mode.eval(at, 0.0, E, X), 0.0, 1e-12);
	EXPECT_NEAR(mode.eval(at, 0.0, E, Y), e_rho, 1e-12);
	EXPECT_NEAR(mode.eval(at, 0.0, E, Z), 0.0, 1e-12);
	EXPECT_NEAR(mode.eval(at, 0.0, H, X), 0.0, 1e-12);
	EXPECT_NEAR(mode.eval(at, 0.0, H, Y), 0.0, 1e-12);
	EXPECT_NEAR(mode.eval(at, 0.0, H, Z), e_rho, 1e-12);
	EXPECT_EQ(113, mode.totalFieldVolume());
	EXPECT_EQ(112, mode.scatteredFieldVolume());
}