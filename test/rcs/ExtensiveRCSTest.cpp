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

#include "TestUtils.h"

namespace maxwell {

using namespace mfem;

class ExtensiveRCSTest : public ::testing::Test {
};

TEST_F(ExtensiveRCSTest, circleTest)
{
	const double dt = 5e-3;
	const int steps = 100000;

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/circle_1m", maxwellCase("2D_RCS"), dt, steps, angles);
}

TEST_F(ExtensiveRCSTest, circleTest_sixmeters)
{
	const double dt = 5e-3;
	const int steps = 100000;

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/circle_salva", maxwellCase("2D_RCS_Salva"), dt, steps, angles);
}

TEST_F(ExtensiveRCSTest, sphereTest)
{
	const double dt = 5e-3;
	const int steps = 100000;

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/sphere", maxwellCase("3D_RCS"), dt, steps, angles);
}
}