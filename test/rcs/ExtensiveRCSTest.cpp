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
#include "components/RCSDataExtractor.h"
#include "TestUtils.h"

namespace maxwell {

using namespace mfem;

class ExtensiveRCSTest : public ::testing::Test {
public:

	template <typename T>
	std::vector<T> linspace(T a, T b, size_t N) {
		T h = (b - a) / static_cast<T>(N-1);
		std::vector<T> xs(N);
		typename std::vector<T>::iterator x;
		T val;
		for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
			*x = val;
		return xs;
	}

	std::vector<SphericalAngles> buildAngleVector(double start_theta, double end_theta, int steps_theta, double start_phi, double end_phi, int steps_phi)
	{
		std::vector<SphericalAngles> res;
		auto theta {linspace(start_theta, end_theta, steps_theta)};
		auto phi {linspace(start_phi, end_phi, steps_phi)};
		for (auto t{0}; t < theta.size(); t++){
			for (auto p {0}; p < phi.size(); p++){
				res.push_back({theta[t], phi[p]});
			}
		}
		return res;
	}


};

TEST_F(ExtensiveRCSTest, circleTest)
{

	std::vector<double> frequencies_manual({ 3e8 / 2.0 / M_PI });

	auto angles{ buildAngleVector(0.0, M_PI, 256, 0.0, M_PI_2, 2) };
	RCSManager rcs("NearToFarFieldExports/circle_1m_O1_Z", maxwellCase("3D_RCS_Sphere_Z"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, circle_1m_O1)
{

	std::vector<double> frequencies_manual({ 3.75e7, 3e8/2.0/M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	// std::vector<double> frequencies_manual({ 3.75e7, 7.5e7 });

	auto angles{ buildAngleVector(M_PI_2, M_PI_2, 1, 0.0, M_PI, 256) };
	RCSManager rcs("NearToFarFieldExports/circle_1m_O1_X", maxwellCase("3D_RCS_Sphere_O1"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, circle_1m_O2)
{

	std::vector<double> frequencies_manual({ 3.75e7, 3e8/2.0/M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	// std::vector<double> frequencies_manual({ 3.75e7 });

	auto angles{ buildAngleVector(0.0, M_PI, 256, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/circle_1m_O2", maxwellCase("3D_RCS_Sphere_O2"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, circle_1m_O3)
{

	std::vector<double> frequencies_manual({ 3.75e7, 3e8/2.0/M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	// std::vector<double> frequencies_manual({ 3.75e7 });

	auto angles{ buildAngleVector(0.0, M_PI, 256, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/circle_1m_O3", maxwellCase("3D_RCS_Sphere_O3"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_1m_O1_Z)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	 std::vector<double> frequencies_manual({ 3e8 });

	auto angles{ buildAngleVector(0.0, M_PI, 256, M_PI, M_PI, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_1m_O1_Z", maxwellCase("3D_RCS_Sphere_1m_O1_Z"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_Box_1m_O1_monostatic)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	//  std::vector<double> frequencies_manual({ 1e8 });
	 auto frequencies_manual = linspace(1e6, 1e9, 200);

	auto angles{ buildAngleVector(M_PI, M_PI, 1, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_Box_1m_G1_O1", maxwellCase("3D_RCS_Sphere_Box_1m_G1_O1"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_Box_1m_O2_monostatic)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	//  std::vector<double> frequencies_manual({ 1e8 });
	auto frequencies_manual = linspace(1e6, 1e9, 200);

	 auto angles{ buildAngleVector(M_PI, M_PI, 1, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_Box_1m_G2_O2", maxwellCase("3D_RCS_Sphere_Box_1m_G2_O2"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_Box_1m_O3_monostatic)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	//  std::vector<double> frequencies_manual({ 1e8 });
	auto frequencies_manual = linspace(1e6, 1e9, 200);

	 auto angles{ buildAngleVector(M_PI, M_PI, 1, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_Box_1m_O3", maxwellCase("3D_RCS_Sphere_Box_1m_O3"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_Box_1m_O1_bistatic)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	std::vector<double> frequencies_manual({ 1.5e8 });

	auto angles{ buildAngleVector(0.0, M_PI, 256, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_Box_1m_G1_O1", maxwellCase("3D_RCS_Sphere_Box_1m_G1_O1"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_Box_1m_O2_bistatic)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	std::vector<double> frequencies_manual({ 1.5e8 });

	auto angles{ buildAngleVector(0.0, M_PI, 256, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_Box_1m_G2_O2", maxwellCase("3D_RCS_Sphere_Box_1m_G2_O2"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, 3D_RCS_Sphere_Box_1m_O3_bistatic)
{

	//std::vector<double> frequencies_manual({ 3.75e7, 3e8 / 2.0 / M_PI, 7.5e7, 1.5e8, 3e8, 6e8, 1e9 });
	std::vector<double> frequencies_manual({ 1e8 });

	auto angles{ buildAngleVector(M_PI, M_PI, 1, 0.0, 0.0, 1) };
	RCSManager rcs("NearToFarFieldExports/sphere_Box_1m_O3", maxwellCase("3D_RCS_Sphere_Box_1m_O3"), frequencies_manual, angles);
}

TEST_F(ExtensiveRCSTest, RCSDataExporter)
{
	EXPECT_NO_THROW(RCSDataExtractor rcsExt("NearToFarFieldExports/circle_1m_O1_X", "3D_RCS_Sphere_O1"));
}

TEST_F(ExtensiveRCSTest, RCSDataExporter_1m)
{
	EXPECT_NO_THROW(RCSDataExtractor rcsExt("NearToFarFieldExports/sphere_1m_O1", "3D_RCS_Sphere_1m_O1"));
}

}