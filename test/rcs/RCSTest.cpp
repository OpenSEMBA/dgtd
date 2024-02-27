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
std::vector<double> linspace(const double min, const double max, const double stepval) 
{
	int steps = int(std::round((max - min) / stepval));
	std::vector<double> res(steps);
	for (int i = 0; i < steps; i++) {
		res[i] = min + stepval * i;
	}
	res.push_back(max);
	return res;
}

};

TEST_F(RCSTest, circleTest)
{
	const double f_min = 1e6;
	const double f_max = 1e9;
	const double f_step = 1e6;

	auto frequency{ linspace(f_min, f_max, f_step) };

	std::vector<double> frequencies_manual({ 30e6, 70e6, 100e6, 200e6, 300e6, 400e6, 500e6, 600e6, 700e6, 800e6, 900e6 });

	std::vector<std::pair<Rho, Phi>> angles{ {0.0, M_PI_2}  }; 

	RCSManager rcs("NearToFarFieldExports/circle", "2D_NTFF_Circle", frequency, angles);
}

}