#include "gtest/gtest.h"

#include "maxwell/Sources.h"

using namespace maxwell;
using namespace mfem;

class TestSources : public ::testing::Test {
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

TEST_F(TestSources, negSpreadInput)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X, -2.0,  1.0, Vector({0.5})));
}

TEST_F(TestSources, negCoeffInput)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X,  2.0, -1.0, Vector({0.5})));
}

TEST_F(TestSources, emptyDevInput)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X,  2.0,  1.0, Vector({   })));
}

TEST_F(TestSources, outOfBoundsDevInput1D)
{
	ASSERT_ANY_THROW(Source(buildOneDimOneMatModel(), E, X,  2.0,  1.0, Vector({ 20.0 })));
}

TEST_F(TestSources, outOfBoundsDevInput2D)
{
	Model model = Model(
		Mesh::MakeCartesian2D(5, 5, Element::Type::QUADRILATERAL), 
		AttributeToMaterial(), 
		AttributeToBoundary());
	ASSERT_ANY_THROW(Source(model, E, X, 2.0, 1.0, Vector({ 20.0, 0.5 })));
	ASSERT_ANY_THROW(Source(model, E, X, 2.0, 1.0, Vector({ 0.5, 20.0 })));
}


