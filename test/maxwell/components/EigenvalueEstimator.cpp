#include <gtest/gtest.h>

#include "TestUtils.h"
#include "components/EigenvalueEstimator.h"
#include "evolution/EvolutionMethods.h"

using namespace maxwell;
using namespace mfem;

using Solver = maxwell::Solver;

class EigenvalueEstimatorTest : public ::testing::Test {
protected:
	Mesh mTri { Mesh::MakeCartesian2D(1, 1,    Element::TRIANGLE     , true) };
	Mesh mQuad{ Mesh::MakeCartesian2D(1, 1,    Element::QUADRILATERAL, true) };
	Mesh mHexa{ Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON) };
	DG_FECollection fec2D{ DG_FECollection(1, 2,BasisType::GaussLegendre) };
	DG_FECollection fec3D{ DG_FECollection(1, 3,BasisType::GaussLegendre) };
	FiniteElementSpace fesTri{  FiniteElementSpace(&mTri,  &fec2D) };
	FiniteElementSpace fesQuad{ FiniteElementSpace(&mQuad, &fec2D) };
	FiniteElementSpace fesHexa{ FiniteElementSpace(&mHexa, &fec3D) };

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
	FiniteElementSpace smFesTri{ FiniteElementSpace(&smTri, &fec2D) };

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
		fesTri,
		Model{
			mTri,
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
		fesTri,
		Model{
			mTri,
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

TEST_F(EigenvalueEstimatorTest, comparePECandSMAconditions)
{
	auto evPEC{EigenvalueEstimator(fesTri, Model{ mTri, AttributeToMaterial{}, AttributeToBoundary{}, AttributeToInteriorConditions{} }, EvolutionOptions{}) };
	auto att_to_bdr{ AttributeToBoundary{{1,BdrCond::SMA}} };
	auto evSMA{ EigenvalueEstimator(fesTri, Model{ mTri, AttributeToMaterial{}, att_to_bdr, AttributeToInteriorConditions{} }, EvolutionOptions{}) };

	auto eigenvalsPEC{ evPEC.getElementMatrix().eigenvalues() };
	auto eigenvalsSMA{ evSMA.getElementMatrix().eigenvalues() };

	std::cout << "PEC" << std::endl;
	std::cout << evPEC.getElementMatrix().eigenvalues() << std::endl;
	std::cout << "SMA" << std::endl;
	std::cout << evSMA.getElementMatrix().eigenvalues() << std::endl;
}