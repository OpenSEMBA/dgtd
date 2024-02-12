#include <gtest/gtest.h>

#include <fftw3.h>

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/RCSExporter.h>
#include <components/RCSManager.h>


namespace maxwell {

using namespace mfem;

class RCSTest : public ::testing::Test {
public:
std::vector<double> linspace(const double min, const double max, const int stepval) 
{
	int steps = int(max / stepval);
	std::vector<double> res(steps);
	for (int i = 0; i < steps; i++) {
		res[i] = min + stepval * i;
	}
	return res;
}

};

TEST_F(RCSTest, circleTest)
{
	const double f_min = 1e6;
	const double f_max = 1e9;
	const int f_step = 1e6;

	auto frequency{ linspace(f_min, f_max, f_step) };

	std::vector<std::pair<Rho, Phi>> angles{ {0.0, 0.0}, {M_PI, 0.0} };
	RCSManager rcs("NearToFarFieldExports/circle", frequency, angles);
}

}