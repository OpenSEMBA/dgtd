#include <gtest/gtest.h>

#include "math/Calculus.h"
#include "components/Types.h"

using namespace maxwell;

class CalculusTest : public ::testing::Test {
};

TEST_F(CalculusTest, crossProduct_basisVectors)
{
	auto i = unitVec(X);
	auto j = unitVec(Y);
	auto k = unitVec(Z);

	// i x j = k
	auto ixj = crossProduct(i, j);
	EXPECT_DOUBLE_EQ(0.0, ixj[0]);
	EXPECT_DOUBLE_EQ(0.0, ixj[1]);
	EXPECT_DOUBLE_EQ(1.0, ixj[2]);

	// j x k = i
	auto jxk = crossProduct(j, k);
	EXPECT_DOUBLE_EQ(1.0, jxk[0]);
	EXPECT_DOUBLE_EQ(0.0, jxk[1]);
	EXPECT_DOUBLE_EQ(0.0, jxk[2]);

	// k x i = j
	auto kxi = crossProduct(k, i);
	EXPECT_DOUBLE_EQ(0.0, kxi[0]);
	EXPECT_DOUBLE_EQ(1.0, kxi[1]);
	EXPECT_DOUBLE_EQ(0.0, kxi[2]);
}

TEST_F(CalculusTest, crossProduct_antiSymmetry)
{
	mfem::Vector a({1.0, 2.0, 3.0});
	mfem::Vector b({4.0, -1.0, 2.0});

	auto axb = crossProduct(a, b);
	auto bxa = crossProduct(b, a);

	for (int i = 0; i < 3; i++) {
		EXPECT_DOUBLE_EQ(-axb[i], bxa[i]);
	}
}

TEST_F(CalculusTest, crossProduct_parallelVectorsGiveZero)
{
	mfem::Vector a({2.0, 3.0, 1.0});
	mfem::Vector b({4.0, 6.0, 2.0}); // 2*a

	auto result = crossProduct(a, b);

	for (int i = 0; i < 3; i++) {
		EXPECT_NEAR(0.0, result[i], 1e-15);
	}
}

TEST_F(CalculusTest, crossProduct_knownResult)
{
	// (1,2,3) x (4,-1,2) = (2*2 - 3*(-1), 3*4 - 1*2, 1*(-1) - 2*4)
	//                     = (4+3, 12-2, -1-8) = (7, 10, -9)
	mfem::Vector a({1.0, 2.0, 3.0});
	mfem::Vector b({4.0, -1.0, 2.0});

	auto r = crossProduct(a, b);

	EXPECT_DOUBLE_EQ( 7.0, r[0]);
	EXPECT_DOUBLE_EQ(10.0, r[1]);
	EXPECT_DOUBLE_EQ(-9.0, r[2]);
}

TEST_F(CalculusTest, unitVec_values)
{
	for (int d = 0; d < 3; d++) {
		auto v = unitVec(d);
		EXPECT_EQ(3, v.Size());
		for (int i = 0; i < 3; i++) {
			if (i == d) {
				EXPECT_DOUBLE_EQ(1.0, v[i]);
			} else {
				EXPECT_DOUBLE_EQ(0.0, v[i]);
			}
		}
	}
}

TEST_F(CalculusTest, minusUnitVec_values)
{
	for (int d = 0; d < 3; d++) {
		auto v = minusUnitVec(d);
		EXPECT_EQ(3, v.Size());
		for (int i = 0; i < 3; i++) {
			if (i == d) {
				EXPECT_DOUBLE_EQ(-1.0, v[i]);
			} else {
				EXPECT_DOUBLE_EQ(0.0, v[i]);
			}
		}
	}
}

TEST_F(CalculusTest, unitVec_orthogonal)
{
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			auto u1 = unitVec(d1);
			auto u2 = unitVec(d2);
			double dot = u1 * u2; // mfem::Vector dot product
			if (d1 == d2) {
				EXPECT_DOUBLE_EQ(1.0, dot);
			} else {
				EXPECT_DOUBLE_EQ(0.0, dot);
			}
		}
	}
}
