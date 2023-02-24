#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>

#include "gtest/gtest.h"

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

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
