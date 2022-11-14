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
TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceNoDir_O1)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(1, elements);
		std::cout << buildMaxwellDGTrace1DEigen(*fes_, {}, 1.0) << std::endl;
		std::cout << buildExpectedJumpDenseMatrix1D(1, elements) << std::endl;
		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {}, 1.0)
			.isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceNoDir_O2)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(2, elements);
		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {}, 1.0)
			.isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceNoDir_O3)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(3, elements);
		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {}, 1.0)
			.isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceNoDir_O4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {}, 1.0)
			.isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXDirO1)
{
	for (int elements = 2; elements < 5; elements++) {
		setFES(1, elements);
		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X}, 1.0)
			.isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));
	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceYDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceYDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceYDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceYDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceZDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceZDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceZDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceZDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXXDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(1, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXXDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(2, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXXDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(3, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXXDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, X}, 1.0).isApprox(buildExpectedJumpDenseMatrix1D(4, elements)));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXYDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXYDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXYDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXYDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Y}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXZDirO1)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(1, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXZDirO2)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(2, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXZDirO3)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(3, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Z}, 1.0).isZero(tol_));

	}
}

TEST_F(MFEMHesthaven1DTrace, MaxwellDGTraceXZDirO4)
{

	for (int elements = 2; elements < 5; elements++) {

		setFES(4, elements);

		EXPECT_TRUE(buildMaxwellDGTrace1DEigen(*fes_, {X, Z}, 1.0).isZero(tol_));

	}
}