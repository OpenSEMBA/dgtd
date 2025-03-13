#include <gtest/gtest.h>

#include <vector>
#include <fftw3.h>

#include <mfem.hpp>

#include <iostream>
#include <filesystem>
#include <math.h>
#include <complex>

#include "components/Types.h"
#include "components/RCSManager.h"
#include "TestUtils.h"

namespace maxwell {

using namespace mfem;

class ExtensiveRCSTest : public ::testing::Test {
};

TEST_F(ExtensiveRCSTest, circleTest)
{

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/circle_1m", maxwellCase("2D_RCS"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, circleTest_sixmeters)
{

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/circle_salva", maxwellCase("2D_RCS_Salva"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, sphereTest)
{

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Phi, Theta>> angles{ {0.0, M_PI_2} };

	RCSManager rcs("NearToFarFieldExports/sphere", maxwellCase("3D_RCS"), frequencies_manual, angles);
}
}