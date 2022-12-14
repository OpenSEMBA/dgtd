#include <gtest/gtest.h>

#include "TestMfemHesthavenFunctions.h"
#include "GlobalFunctions.h"

using namespace maxwell;
using namespace mfem;

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

	static FluxCoefficient buildAverageOnly() { return { 1.0 }; }
	static FluxCoefficient buildJumpOnly() { return { 1.0 }; }

};

//TEST_F(MFEMHesthaven1DTrace, StrongFluxOperator_O1)
//{
//	auto mInv = buildInverseMassMatrixEigen(*fes_);
//	auto flux = buildNormalPECFluxOperator1D(*fes_, {X});
//
//}

TEST_F(MFEMHesthaven1DTrace, DGTraceAverageOnly_O1)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(1, elements);
		EXPECT_TRUE(buildDGTraceAverage1DEigen(*fes_, buildAverageOnly())
			.isApprox(buildExpectedAverageDenseMatrix1D(1, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, DGTraceAverageOnly_O2)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(2, elements);
		EXPECT_TRUE(buildDGTraceAverage1DEigen(*fes_, buildAverageOnly())
			.isApprox(buildExpectedAverageDenseMatrix1D(2, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, DGTraceAverageOnly_O3)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(3, elements);
		EXPECT_TRUE(buildDGTraceAverage1DEigen(*fes_, buildAverageOnly())
			.isApprox(buildExpectedAverageDenseMatrix1D(3, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, DGTraceAverageOnly_O4)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(4, elements);
		EXPECT_TRUE(buildDGTraceAverage1DEigen(*fes_, buildAverageOnly())
			.isApprox(buildExpectedAverageDenseMatrix1D(4, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, DGTraceJumpOnly_O1)
{

	for (int elements = 2; elements < 5; elements++) {
		setFES(1, elements);
		EXPECT_TRUE(buildDGTraceJump1DEigen(*fes_, buildJumpOnly())
			.isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, DGTraceJumpOnly_O2)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(2, elements);

		EXPECT_TRUE(
			buildDGTraceJump1DEigen(*fes_, buildJumpOnly())
			.isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, DGTraceJumpOnly_O3)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(3, elements);
		EXPECT_TRUE(buildDGTraceJump1DEigen(*fes_, buildJumpOnly())
			.isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));
	}
}
TEST_F(MFEMHesthaven1DTrace, DGTraceJumpOnly_O4)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(4, elements);
		EXPECT_TRUE(buildDGTraceJump1DEigen(*fes_, buildJumpOnly())
			.isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));
	}
}