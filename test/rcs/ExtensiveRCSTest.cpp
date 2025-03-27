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
public:
	std::vector<SphericalAngles> buildAngleVector(double start_phi, double end_phi, int steps_phi, double start_theta, double end_theta, int steps_theta)
	{
		std::vector<SphericalAngles> res;
		auto phi_incr{ (end_phi - start_phi) / steps_phi };
		auto the_incr{ (end_theta - start_theta) / steps_theta };
		for (auto p{ 0 }; p < steps_phi; p++) {
			for (auto t{ 0 }; t < steps_theta; t++) {
				res.push_back({ start_phi + p * phi_incr, start_theta + t * the_incr });
			}
		}
		return res;
	}

};

TEST_F(ExtensiveRCSTest, circleTest)
{

	std::vector<double> frequencies_manual({ 3e8 });

	auto angles{ buildAngleVector(M_PI, 2.0 * M_PI, 128, M_PI_2, M_PI_2, 1) };
	RCSManager rcs("NearToFarFieldExports/circle_1m_O1", maxwellCase("3D_RCS_Sphere"), frequencies_manual, angles);
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
