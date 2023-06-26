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
	Mesh mQuad{ Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL, true) };
	DG_FECollection fec{ DG_FECollection(1, 2,BasisType::GaussLegendre) };
	FiniteElementSpace fesTri{ FiniteElementSpace(&mTri, &fec) };
	FiniteElementSpace fesQuad{ FiniteElementSpace(&mQuad, &fec) };

};

TEST_F(EigenvalueEstimatorTest, verifyNoThrow_Quad)
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

TEST_F(EigenvalueEstimatorTest, verifyNoThrow_SubMesh_Tri)
{
	mTri.SetAttribute(0, 300);
	Array<int> atts(1); atts[0] = 300;
	auto smTri(SubMesh::CreateFromDomain(mTri, atts));
	smTri.SetAttribute(0, 1);
	FiniteElementSpace smFesTri{ FiniteElementSpace(&smTri, &fec) };

	EXPECT_NO_THROW(EigenvalueEstimator(
		smFesTri, 
		Model(
			smTri, 
			AttributeToMaterial{},
			AttributeToBoundary{},
			AttributeToInteriorConditions{}),
		EvolutionOptions{})
	);
}

TEST_F(EigenvalueEstimatorTest, printMatrix_PEC)
{
	EigenvalueEstimator ev(
		fesQuad,
		Model{
			mQuad,
			AttributeToMaterial{},
			AttributeToBoundary{},
			AttributeToInteriorConditions{} },
		EvolutionOptions{});

	auto mat{ ev.getElementMatrix() };
	
	ASSERT_GE(mat.rows(), 0);
	ASSERT_GE(mat.cols(), 0);

	std::cout << mat << std::endl;
	std::cout << std::flush;

}

TEST_F(EigenvalueEstimatorTest, printMatrix_SMA)
{

	EigenvalueEstimator ev(
		fesQuad,
		Model{
			mQuad,
			AttributeToMaterial{},
			AttributeToBoundary{ {1,BdrCond::SMA} },
			AttributeToInteriorConditions{} },
			EvolutionOptions{});

	auto mat{ ev.getElementMatrix() };

	ASSERT_GE(mat.rows(), 0);
	ASSERT_GE(mat.cols(), 0);

	std::cout << mat << std::endl;
	std::cout << std::flush;

}