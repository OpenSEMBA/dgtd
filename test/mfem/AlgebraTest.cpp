#include <gtest/gtest.h>

#include <iostream>
#include <mfem.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <math.h>

using namespace mfem;
class AlgebraTest : public ::testing::Test 
{
public:

	double speed_of_light{ 299792458.0 };
	double func_exp_real_part_2D(const Vector& x, const double freq, const double phi)
	{
		//angulo viene dado por x[0], x[1] y 0.0, 0.0. No es el angulo donde observo, es el angulo que forma el punto y el angulo de observacion en un sistema centrado en el punto.
		auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
		return cos(2.0 * M_PI * (speed_of_light / freq) * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
	}

	double func_exp_imag_part_2D(const Vector& x, const double freq, const double phi)
	{
		auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
		return sin(2.0 * M_PI * (speed_of_light / freq) * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
	}

	std::complex<double> func_exp_2D(const Vector& x, const double freq, const double phi)
	{
		auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
		return exp(std::complex(0.0, 2.0 * M_PI * (speed_of_light / freq) * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle)));
	}
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

TEST_F(AlgebraTest, custom_rcs_exponential_function)
{
	Vector vec({ 1.5, 1.5 });
	double tol{ 1e-5 };

	EXPECT_NEAR(func_exp_real_part_2D(vec, 100000.0, 0.0),  0.757829, tol);
	EXPECT_NEAR(func_exp_imag_part_2D(vec, 100000.0, 0.0), -0.652453, tol);
	EXPECT_NEAR(func_exp_real_part_2D(vec, 100000.0, M_PI), 0.757829, tol);
	EXPECT_NEAR(func_exp_imag_part_2D(vec, 100000.0, M_PI), 0.652453, tol);

	EXPECT_NEAR(func_exp_2D(vec, 100000.0, 0.0).real(), func_exp_real_part_2D(vec, 100000.0, 0.0), tol);
	EXPECT_NEAR(func_exp_2D(vec, 100000.0, 0.0).imag(), func_exp_imag_part_2D(vec, 100000.0, 0.0), tol);
	EXPECT_NEAR(func_exp_2D(vec, 100000.0, M_PI).real(), func_exp_real_part_2D(vec, 100000.0, M_PI), tol);
	EXPECT_NEAR(func_exp_2D(vec, 100000.0, M_PI).imag(), func_exp_imag_part_2D(vec, 100000.0, M_PI), tol);
}