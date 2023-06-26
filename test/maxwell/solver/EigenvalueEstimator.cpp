#include <gtest/gtest.h>

#include "TestUtils.h"
#include "evolution/EigenvalueEstimator.h"
#include "evolution/EvolutionMethods.h"

using namespace maxwell;
using namespace mfem;

using Solver = maxwell::Solver;

class EigenvalueEstimatorTest : public ::testing::Test {
protected:
	
	Mesh mTri{ Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE, true) };
	Mesh mQuad{ Mesh::MakeCartesian2D(1,1, Element::QUADRILATERAL, true) };
	DG_FECollection fec{ DG_FECollection(2,1,BasisType::GaussLegendre) };
	FiniteElementSpace fesTri{ FiniteElementSpace(&mTri, &fec) };
	FiniteElementSpace fesQuad{ FiniteElementSpace(&mQuad, &fec) };
};

TEST_F(EigenvalueEstimatorTest, assemblerNoThrow_Quad)
{
	EXPECT_NO_THROW(EigenvalueEstimator(
		fesQuad, 
		Model{ 
			mQuad, 
			AttributeToMaterial{}, 
			AttributeToBoundary{}, 
			AttributeToInteriorConditions{} }, 
		EvolutionOptions{})
	);
}