#include "gtest/gtest.h"

#include "maxwell/Sources.h"

using namespace maxwell;

class TestMaxwellSources : public ::testing::Test {
protected:

	AttributeToBoundary buildAttrToBdrMap1D(const BdrCond& bdrL, const BdrCond& bdrR)
	{
		return {
			{1, bdrL},
			{2, bdrR}
		};
	}

	Model buildOneDimOneMatModel(
		const int meshIntervals = 51,
		const BdrCond& bdrL = BdrCond::PEC,
		const BdrCond& bdrR = BdrCond::PEC) {

		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), AttributeToMaterial(), buildAttrToBdrMap1D(bdrL, bdrR));
	}

};

TEST_F(TestMaxwellSources, negSpreadInput)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X, -2.0,  1.0, Vector({0.5})));
}

TEST_F(TestMaxwellSources, negCoeffInput)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X,  2.0, -1.0, Vector({0.5})));
}

TEST_F(TestMaxwellSources, emptyDevInput)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X,  2.0,  1.0, Vector({   })));
}

TEST_F(TestMaxwellSources, outOfBoundsDevInput1D)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X,  2.0,  1.0, Vector({ 20.0 })));
}

TEST_F(TestMaxwellSources, outOfBoundsDevInput2D)
{
	Model model = Model(
		Mesh::MakeCartesian2D(5, 5, Element::Type::QUADRILATERAL), 
		AttributeToMaterial(), 
		AttributeToBoundary());
	ASSERT_ANY_THROW(Source(model, E, X, 2.0, 1.0, Vector({ 20.0, 0.5 })));
	ASSERT_ANY_THROW(Source(model, E, X, 2.0, 1.0, Vector({ 0.5, 20.0 })));
}


