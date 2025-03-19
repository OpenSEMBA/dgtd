#include "math/Function.h"
#include <gtest/gtest.h>

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