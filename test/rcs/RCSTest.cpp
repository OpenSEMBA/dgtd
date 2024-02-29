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

namespace maxwell {

using namespace mfem;

class RCSTest : public ::testing::Test {
public:
std::vector<double> linspace(const double max, const int steps) 
{
	std::vector<double> res(steps+1);
	const double step = double(max / steps);
	for (int i = 0; i <= steps; i++) {
		res[i] = step * i;
	}
	return res;
}

};

TEST_F(RCSTest, circleTest)
{
	const double dt = 5e-3;
	const double f_max = 2.0 / dt;
	const int steps = 100000;

	auto frequency{ linspace(f_max, steps) };

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Rho, Phi>> angles{ {0.0, M_PI_2}  }; 

	RCSManager rcs("NearToFarFieldExports/circle", "2D_NTFF_Circle", frequency, angles);
}

}