#include <gtest/gtest.h>

#include "mfemExtension/IntegratorFunctions.h"
#include "components/Types.h"

using namespace maxwell;
using namespace maxwell::mfemExtension;
using namespace mfem;

class IntegratorFunctionsTest : public ::testing::Test {
};

// --- buildNormalTerm ---

TEST_F(IntegratorFunctionsTest, buildNormalTerm_1D)
{
	Vector nor({1.0});
	EXPECT_DOUBLE_EQ( 1.0, buildNormalTerm(nor, X));
	EXPECT_DOUBLE_EQ( 0.0, buildNormalTerm(nor, Y)); // padded to zero
	EXPECT_DOUBLE_EQ( 0.0, buildNormalTerm(nor, Z)); // padded to zero
}

TEST_F(IntegratorFunctionsTest, buildNormalTerm_2D)
{
	Vector nor({3.0, -2.0});
	EXPECT_DOUBLE_EQ( 3.0, buildNormalTerm(nor, X));
	EXPECT_DOUBLE_EQ(-2.0, buildNormalTerm(nor, Y));
	EXPECT_DOUBLE_EQ( 0.0, buildNormalTerm(nor, Z));
}

TEST_F(IntegratorFunctionsTest, buildNormalTerm_3D)
{
	Vector nor({1.0, 2.0, 3.0});
	EXPECT_DOUBLE_EQ(1.0, buildNormalTerm(nor, X));
	EXPECT_DOUBLE_EQ(2.0, buildNormalTerm(nor, Y));
	EXPECT_DOUBLE_EQ(3.0, buildNormalTerm(nor, Z));
}

TEST_F(IntegratorFunctionsTest, buildNormalTerm_negativeComponents)
{
	Vector nor({-0.5, 0.0, 0.7});
	EXPECT_DOUBLE_EQ(-0.5, buildNormalTerm(nor, X));
	EXPECT_DOUBLE_EQ( 0.0, buildNormalTerm(nor, Y));
	EXPECT_DOUBLE_EQ( 0.7, buildNormalTerm(nor, Z));
}

// --- buildFaceMatrix ---

TEST_F(IntegratorFunctionsTest, buildFaceMatrix_diagonalBlock_adds)
{
	// offsetRow == offsetCol => += w * shapeA(i) * shapeB(j)
	DenseMatrix elmat(4, 4);
	elmat = 0.0;

	Vector shapeA(2);  shapeA  = 1.0;
	Vector shapeB(2);  shapeB  = 1.0;

	buildFaceMatrix(2.0, 2, 2, 0, 0, shapeA, shapeB, elmat);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			EXPECT_DOUBLE_EQ(2.0, elmat(i, j));
		}
	}
	// Off-diagonal block should be untouched
	for (int i = 0; i < 2; i++) {
		for (int j = 2; j < 4; j++) {
			EXPECT_DOUBLE_EQ(0.0, elmat(i, j));
		}
	}
}

TEST_F(IntegratorFunctionsTest, buildFaceMatrix_offDiagonalBlock_subtracts)
{
	// offsetRow != offsetCol => -= w * shapeA(i) * shapeB(j)
	DenseMatrix elmat(4, 4);
	elmat = 0.0;

	Vector shapeA(2);  shapeA  = 1.0;
	Vector shapeB(2);  shapeB  = 1.0;

	buildFaceMatrix(3.0, 2, 2, 0, 2, shapeA, shapeB, elmat);

	for (int i = 0; i < 2; i++) {
		for (int j = 2; j < 4; j++) {
			EXPECT_DOUBLE_EQ(-3.0, elmat(i, j));
		}
	}
	// The (0,0) block should be untouched
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			EXPECT_DOUBLE_EQ(0.0, elmat(i, j));
		}
	}
}

TEST_F(IntegratorFunctionsTest, buildFaceMatrix_accumulatesCorrectly)
{
	// Two calls add on top of each other
	DenseMatrix elmat(2, 2);
	elmat = 0.0;

	Vector shapeA(2);  shapeA  = 1.0;
	Vector shapeB(2);  shapeB  = 1.0;

	buildFaceMatrix(1.0, 2, 2, 0, 0, shapeA, shapeB, elmat);
	buildFaceMatrix(2.0, 2, 2, 0, 0, shapeA, shapeB, elmat);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			EXPECT_DOUBLE_EQ(3.0, elmat(i, j));
		}
	}
}

TEST_F(IntegratorFunctionsTest, buildFaceMatrix_nonUniformShapes)
{
	DenseMatrix elmat(2, 2);
	elmat = 0.0;

	Vector shapeA(2);
	shapeA(0) = 2.0;
	shapeA(1) = 3.0;
	Vector shapeB(2);
	shapeB(0) = 1.0;
	shapeB(1) = 4.0;

	buildFaceMatrix(1.0, 2, 2, 0, 0, shapeA, shapeB, elmat);

	EXPECT_DOUBLE_EQ(2.0 * 1.0, elmat(0, 0));
	EXPECT_DOUBLE_EQ(2.0 * 4.0, elmat(0, 1));
	EXPECT_DOUBLE_EQ(3.0 * 1.0, elmat(1, 0));
	EXPECT_DOUBLE_EQ(3.0 * 4.0, elmat(1, 1));
}

// --- calculateBetaTerm ---

TEST_F(IntegratorFunctionsTest, calculateBetaTerm_noDirections_returnsAlpha)
{
	Vector nor({0.0, 1.0, 0.0});
	std::vector<Direction> dirs;
	EXPECT_DOUBLE_EQ(5.0, calculateBetaTerm(nor, dirs, 5.0));
}

TEST_F(IntegratorFunctionsTest, calculateBetaTerm_oneDirection)
{
	Vector nor({0.0, 0.8, 0.0});
	std::vector<Direction> dirs = {Y};
	EXPECT_DOUBLE_EQ(3.0 * 0.8, calculateBetaTerm(nor, dirs, 3.0));
}

TEST_F(IntegratorFunctionsTest, calculateBetaTerm_twoDirections)
{
	Vector nor({0.6, 0.8, 0.0});

	std::vector<Direction> dirsXX = {X, X};
	EXPECT_NEAR(1.0 * 0.6 * 0.6, calculateBetaTerm(nor, dirsXX, 1.0), 1e-15);

	std::vector<Direction> dirsXY = {X, Y};
	EXPECT_NEAR(1.0 * 0.6 * 0.8, calculateBetaTerm(nor, dirsXY, 1.0), 1e-15);

	std::vector<Direction> dirsYY = {Y, Y};
	EXPECT_NEAR(1.0 * 0.8 * 0.8, calculateBetaTerm(nor, dirsYY, 1.0), 1e-15);
}

TEST_F(IntegratorFunctionsTest, calculateBetaTerm_threeDirections_throws)
{
	Vector nor({1.0, 0.0, 0.0});
	std::vector<Direction> dirs = {X, X, X};
	EXPECT_ANY_THROW(calculateBetaTerm(nor, dirs, 1.0));
}

TEST_F(IntegratorFunctionsTest, calculateBetaTerm_zeroNormal)
{
	Vector nor({0.0, 0.0, 0.0});
	std::vector<Direction> dirs = {X};
	EXPECT_DOUBLE_EQ(0.0, calculateBetaTerm(nor, dirs, 2.0));
}

// --- setNormalVector1D ---

TEST_F(IntegratorFunctionsTest, setNormalVector1D_leftBoundary)
{
	// x = 0 => 2*0 - 1 = -1 (outward normal at left boundary)
	IntegrationPoint ip;
	ip.x = 0.0;

	auto nor = setNormalVector1D(1, ip);
	ASSERT_EQ(1, nor.Size());
	EXPECT_DOUBLE_EQ(-1.0, nor[0]);
}

TEST_F(IntegratorFunctionsTest, setNormalVector1D_rightBoundary)
{
	// x = 1 => 2*1 - 1 = 1 (outward normal at right boundary)
	IntegrationPoint ip;
	ip.x = 1.0;

	auto nor = setNormalVector1D(1, ip);
	ASSERT_EQ(1, nor.Size());
	EXPECT_DOUBLE_EQ(1.0, nor[0]);
}

TEST_F(IntegratorFunctionsTest, setNormalVector1D_midPoint)
{
	// x = 0.5 => 2*0.5 - 1 = 0
	IntegrationPoint ip;
	ip.x = 0.5;

	auto nor = setNormalVector1D(1, ip);
	ASSERT_EQ(1, nor.Size());
	EXPECT_DOUBLE_EQ(0.0, nor[0]);
}
