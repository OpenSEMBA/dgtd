#include "gtest/gtest.h"

#include <iostream>
#include <mfem.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <math.h>

using namespace mfem;
class AlgebraTest : public ::testing::Test {

};

std::complex<double> operator"" _i(long double x)
{
	return std::complex<double>(0.0, x);
}

TEST_F(AlgebraTest, calcRealEigenvalues)
{

	Eigen::Matrix3d matrix{
		{1.0, 0.0, 0.0},
		{0.0, 2.0, 0.0},
		{0.0, 0.0,-3.0}
	};

	Eigen::Vector<double, 3> expectedEVs{
		{1.0, 2.0, -3.0}
	};

	EXPECT_EQ(matrix.eigenvalues(), expectedEVs);

}

TEST_F(AlgebraTest, calcComplexEigenvalues)
{

	Eigen::Matrix3cd matrix{
	{1.0 + 1.0_i, 0.0        , 0.0        },
	{0.0        , 2.0 + 2.0_i, 0.0        },
	{0.0        , 0.0        ,-3.0 + 3.0_i}
	};

	Eigen::Vector<std::complex<double>, 3> expectedEVs{
	{1.0 + 1.0_i, 2.0 + 2.0_i, -3.0 + 3.0_i}
	};

	EXPECT_EQ(matrix.eigenvalues(), expectedEVs);

}