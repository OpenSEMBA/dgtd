#include "gtest/gtest.h"

#include "maxwell/Sources.h"

using namespace maxwell;
using namespace mfem;

using Position = GaussianInitialField::Position;

class TestSources : public ::testing::Test {
protected:

};

TEST_F(TestSources, negSpreadInput)
{
	ASSERT_ANY_THROW(GaussianInitialField(E, X, -2.0,  1.0, Position({ 0.5 })));
}

TEST_F(TestSources, negNormalizationInput)
{
	ASSERT_ANY_THROW(GaussianInitialField(E, X,  2.0, -1.0, Position({0.5})));
}


