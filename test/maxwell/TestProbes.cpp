#include "gtest/gtest.h"

#include "maxwell/Probes.h"

using namespace maxwell;
using mfem::DenseMatrix;

class TestProbes : public ::testing::Test {
};

TEST_F(TestProbes, integPointDiffDimsVector)
{
	auto pointVec = std::vector<std::vector<double>>({ {0.0, 0.5}, {0.5}, {1.0, 0.5, 1.0} });
	ASSERT_ANY_THROW( PointsProbe(E, X, pointVec) );
}

TEST_F(TestProbes, integPointEmptySubvectors)
{
	auto pointVec = std::vector<std::vector<double>>({ {},{} });
	ASSERT_ANY_THROW(PointsProbe(E, X, pointVec));
}