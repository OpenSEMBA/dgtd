#include <gtest/gtest.h>

#include <iostream>
#include <mfem.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

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
	Eigen::Matrix2d matrix{
		{  3, -2},
		{  4, -1}
	};

	Eigen::Vector<std::complex<double>, 2> expectedEVs{
		{1.0 + 2.0_i, 1.0 - 2.0_i}
	};

	EXPECT_EQ(matrix.eigenvalues(), expectedEVs);

}

TEST_F(AlgebraTest, checkSparseMatrixVectorProduct)
{
	int a{ 4 }, b{ 4 };
	Eigen::SparseMatrix<double> sparse(a,b);
	sparse.insert(0, 0) = 1.0;
	sparse.insert(1, 0) = 2.0;
	sparse.insert(0, 2) = 3.0;
	sparse.insert(2, 3) = 4.0;
	sparse.makeCompressed();
	int c = 0;
	Eigen::Vector4d vec({ 1.0,1.0,1.0,1.0 });
	
	EXPECT_NO_THROW( sparse * vec );
}

TEST_F(AlgebraTest, checkEigenMatrixVectorProduct)
{
	int a{ 2 }, b{ 2 };
	Eigen::SparseMatrix<double> sparse(a, b);
	sparse.insert(0, 0) = 2;
	sparse.insert(0, 1) = -12;
	sparse.insert(1, 0) = 1;
	sparse.insert(1, 1) = -5;
	Eigen::VectorXd it(2);
	it[0] = 1;
	it[1] = 1;

	EXPECT_NO_THROW( sparse * it);
}