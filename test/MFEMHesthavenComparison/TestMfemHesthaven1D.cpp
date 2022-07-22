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
		mesh_ = Mesh::MakeCartesian1D(1);
		fec_ = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
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

	EXPECT_TRUE(buildMassMatrixEigen(fes_).isApprox(expMat, tol_));
}

TEST_F(MFEMHesthaven1D, checkInverseMassMatrix)
{
	
	auto expMat = Eigen::Matrix2d{
		{4.0, -2.0},
		{-2.0, 4.0}, };

	EXPECT_TRUE(buildInverseMassMatrixEigen(fes_).isApprox(expMat, tol_));

}

TEST_F(MFEMHesthaven1D, checkStiffnessMatrixO1)
{
	
	auto expMat = Eigen::Matrix2d{
			{-0.5, 0.5},
			{-0.5, 0.5} };

	EXPECT_TRUE(buildStiffnessMatrixEigen(fes_).isApprox(expMat, tol_));

}
TEST_F(MFEMHesthaven1D, checkDOperatorO1)
{

	Eigen::Matrix2d DMatrix = 0.5 * buildInverseMassMatrixEigen(fes_) * buildStiffnessMatrixEigen(fes_);

	std::cout << DMatrix << std::endl;

	auto expMat = Eigen::Matrix2d{
			{-0.5, 0.5},
			{-0.5, 0.5} };

	EXPECT_TRUE(DMatrix.isApprox(expMat, tol_));
}

TEST_F(MFEMHesthaven1D, checkDOperatorO2)
{

	std::unique_ptr<FiniteElementSpace> fes = buildFiniteElementSpace(2);

	Eigen::Matrix3d DMatrix = 0.5 * buildInverseMassMatrixEigen(fes) * buildStiffnessMatrixEigen(fes);
	auto expMat = Eigen::Matrix3d{
			{-1.5, 2.0,-0.5},
			{-0.5, 0.0, 0.5},
			{ 0.5,-2.0, 1.5} };
	EXPECT_TRUE(DMatrix.isApprox(expMat, tol_));
}

TEST_F(MFEMHesthaven1D, checkDOperatorO4)
{

	std::unique_ptr<FiniteElementSpace> fes = buildFiniteElementSpace(4);

	Eigen::Matrix<double,5,5> DMatrix = 0.5 * buildInverseMassMatrixEigen(fes) * buildStiffnessMatrixEigen(fes);
	auto expMat = Eigen::MatrixXd{
			{-5.000000000000000e+00, 6.756502488724238e+00,-2.666666666666666e+00, 1.410164177942427e+00,-5.000000000000000e-01},
			{-1.240990253030983e+00, 0.000000000000000e+00, 1.745743121887938e+00,-7.637626158259730e-01, 2.590097469690174e-01},
			{ 3.750000000000002e-01,-1.336584577695454e+00, 0.000000000000000e+00, 1.336584577695453e+00,-3.750000000000002e-01},
			{-2.590097469690172e-01, 7.637626158259738e-01,-1.745743121887938e+00, 0.000000000000000e+00, 1.240990253030984e+00},
			{ 5.000000000000000e+00,-1.410164177942426e+00, 2.666666666666663e+00,-6.756502488724239e+00, 5.000000000000001e+00} };

	EXPECT_TRUE(DMatrix.isApprox(expMat, tol_));

}


TEST_F(MFEMHesthaven1D, checkStrongFluxOperator)
{

	std::unique_ptr<FiniteElementSpace> fes = buildFiniteElementSpace(2);

	auto mInv = buildInverseMassMatrixEigen(fes);
	std::cout << mInv << std::endl;
	auto flux = buildNormalPECFluxOperator1D(fes, std::vector<maxwell::Direction>{maxwell::Direction::X});
	std::cout << flux << std::endl;

}