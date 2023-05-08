#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>

#include "gtest/gtest.h"

#include <complex>
#include <cmath>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues>

using namespace mfem;

class EigenTest : public ::testing::Test {
};

TEST_F(EigenTest, sparseSavingMethods) 
{
	Eigen::SparseMatrix<double> sparse(5,5);
	sparse.insert(1, 2) = 3.0;
	sparse.insert(2, 4) = 5.0;
	sparse.insert(0, 4) = 1.0;
	Eigen::saveMarket(sparse, "SparseMatrix.mtx");
}

TEST_F(EigenTest, calculatingEigenvalues)
{
	std::complex<double> complex(3.0,2.0);
	std::complex<double> complex2(4.0,3.0);
	Eigen::SparseMatrix<std::complex<double>> sparse(5, 5);
	sparse.insert(0, 0) = 1.0;
	sparse.insert(1, 1) = 2.0;
	sparse.insert(2, 2) = complex;
	sparse.insert(3, 3) = complex2;
	sparse.insert(4, 4) = 5.0;
	Eigen::Vector<std::complex<double>, 5> expected{ 1.0, 2.0, complex, complex2, 5.0 };
	EXPECT_TRUE(expected.isApprox(sparse.toDense().eigenvalues()));

}
