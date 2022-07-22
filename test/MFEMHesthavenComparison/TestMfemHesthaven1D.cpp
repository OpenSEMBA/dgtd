#include "gtest/gtest.h"

#include "mfem.hpp"
#include "TestMfemHesthavenFunctions.cpp"

#include <fstream>
#include <iostream>
#include <../../maxwell/src/maxwell/Types.h>
#include <Eigen/Dense>

using namespace mfem;

class MFEMHesthaven1D : public ::testing::Test {
protected:

	void SetUp() override {
		Mesh mesh_ = Mesh::MakeCartesian1D(1);
		std::unique_ptr<FiniteElementCollection> fec_ = std::make_unique<DG_FECollection>(1, 1);
		std::unique_ptr<FiniteElementSpace> fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};

TEST_F(MFEMHesthaven1D, checkMassMatrix)
{
	auto expMat = Eigen::Matrix2d{
		{2.0 / 6.0, 1.0 / 6.0},
		{1.0 / 6.0, 2.0 / 6.0}, };

	EXPECT_TRUE(buildMassMatrixEigen(fes_.get()).isApprox(expMat, tol_));
}

TEST_F(MFEMHesthaven1D, checkInverseMassMatrix)
{
	auto expMat = Eigen::Matrix2d{
		{4.0, -2.0},
		{-2.0, 4.0}, };

	EXPECT_TRUE(buildInverseMassMatrixEigen(fes_.get()).isApprox(expMat, tol_));

}

TEST_F(MFEMHesthaven1D, checkStiffnessMatrixO1)
{
	auto expMat = Eigen::Matrix2d{
			{-0.5, 0.5},
			{-0.5, 0.5} };

	EXPECT_TRUE(buildStiffnessMatrixEigen(fes_.get()).isApprox(expMat, tol_));

}
TEST_F(MFEMHesthaven1D, checkDOperatorO1)
{

	auto DMatrix = 0.5 * buildInverseMassMatrixEigen(fes_.get()) * buildStiffnessMatrixEigen(fes_.get());
	auto expMat = Eigen::Matrix2d{
			{-0.5, 0.5},
			{-0.5, 0.5} };

	EXPECT_TRUE(DMatrix.isApprox(expMat, tol_));
}

TEST_F(MFEMHesthaven1D, checkDOperatorO2)
{

	std::unique_ptr<FiniteElementSpace> fes = buildFiniteElementSpace(2);

	auto DMatrix = 0.5 * buildInverseMassMatrixEigen(fes_.get()) * buildStiffnessMatrixEigen(fes_.get());
	auto expMat = Eigen::Matrix3d{
			{-1.5, 2.0,-0.5},
			{-0.5, 0.0, 0.5},
			{ 0.5,-2.0, 1.5} };
	EXPECT_TRUE(DMatrix.isApprox(expMat, tol_));
}

TEST_F(MFEMHesthaven1D, checkDOperatorO4)
{

		std::unique_ptr<FiniteElementSpace> fes = buildFiniteElementSpace(4);

		auto MISCalcOld = 0.5 * buildInverseMassMatrixEigen(fes.get()) * buildStiffnessMatrixEigen(fes.get());
		auto expMat = Eigen::MatrixXd{
				{-5.00, 6.76,-2.67, 1.41,-0.50},
				{-1.24, 0.00, 1.75,-0.76, 0.26},
				{ 0.38,-1.34, 0.00, 1.34,-0.38},
				{-0.26, 0.76,-1.75, 0.00, 1.24},
				{ 0.50,-1.41, 2.67,-6.76, 5.00} };

		EXPECT_TRUE(MISCalcOld.isApprox(expMat, 1e-1));

}


TEST_F(MFEMHesthaven1D, checkStrongFluxOperator)
{

	std::unique_ptr<FiniteElementSpace> fes = buildFiniteElementSpace(2);

	auto mInv = buildInverseMassMatrixEigen(fes.get());
	std::cout << mInv << std::endl;
	auto flux = buildNormalPECFluxOperator1D(fes.get(), std::vector<maxwell::Direction>{maxwell::Direction::X});
	std::cout << flux << std::endl;

}