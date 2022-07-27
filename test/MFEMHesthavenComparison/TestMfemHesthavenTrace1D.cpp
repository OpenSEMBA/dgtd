#pragma once

#include "TestMfemHesthavenFunctions.h"
#include "../TestGlobalFunctions.h"

class MFEMHesthaven1DTrace : public ::testing::Test {
protected:

	void SetUp() override
	{
		mesh_ = Mesh::MakeCartesian1D(1);
		fec_ = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES(const int order, const int elements = 1)
	{
		mesh_ = Mesh::MakeCartesian1D(elements);
		fec_ = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};


TEST_F(MFEMHesthaven1DTrace, checkStrongFluxOperator)
{

	setFES(2);

	auto mInv = buildInverseMassMatrixEigen(fes_);
	std::cout << mInv << std::endl;
	auto flux = buildNormalPECFluxOperator1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X});
	std::cout << flux << std::endl;

}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceAverageOnlyMatrixO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(1.0, 0.0)).isApprox(buildExpectedAverageDenseMatrix1D(1, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceAverageOnlyMatrixO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(1.0, 0.0)).isApprox(buildExpectedAverageDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceAverageOnlyMatrixO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(1.0, 0.0)).isApprox(buildExpectedAverageDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceAverageOnlyMatrixO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(1.0, 0.0)).isApprox(buildExpectedAverageDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceJumpOnlyMatrixO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(0.0, 1.0)).isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceJumpOnlyMatrixO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(0.0, 1.0)).isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceJumpOnlyMatrixO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(0.0, 1.0)).isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkDGTraceJumpOnlyMatrixO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenDGTrace1D(fes_.get(), std::make_pair<double, double>(0.0, 1.0)).isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}


TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceNoDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceNoDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceNoDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceNoDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceYDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceYDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceYDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceYDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceZDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceZDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceZDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceZDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXXDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXXDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXXDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXXDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXYDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXYDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXYDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXYDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXZDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXZDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXZDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, checkMaxwellDGTraceXZDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildEigenMaxwellDGTrace1D(fes_.get(), std::vector<maxwell::Direction>{maxwell::Direction::X, maxwell::Direction::Z}, 1.0).isZero(tol_));

	}
}