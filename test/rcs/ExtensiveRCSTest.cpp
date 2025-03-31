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

	std::vector<SphericalAngles> buildAngleVector(double start_phi, double end_phi, int steps_phi, double start_theta, double end_theta, int steps_theta)
	{
		std::vector<SphericalAngles> res;
		auto phi {linspace(start_phi, end_phi, steps_phi)};
		auto theta {linspace(start_theta, end_theta, steps_theta)};
		for (auto p {0}; p < phi.size(); p++){
			for (auto t{0}; t < theta.size(); t++){
				res.push_back({theta[t], phi[p]});
			}
		}
		return res;
	}


};

TEST_F(ExtensiveRCSTest, circleTest)
{

	std::vector<double> frequencies_manual({ 3e8 });

	auto angles{ buildAngleVector(0.0, M_PI_2, 3, 0.0, M_PI, 128) };
	RCSManager rcs("NearToFarFieldExports/circle_1m_O1_Z", maxwellCase("3D_RCS_Sphere_Z"), frequencies_manual, angles);
}

}