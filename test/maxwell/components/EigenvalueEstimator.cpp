#include "components/EigenvalueEstimator.h"

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
	auto map = GeomTagToBoundary{};
	EXPECT_NO_THROW(EigenvalueEstimator(
		fesQuad, 
		Model{ 
			mQuad, 
			GeomTagToMaterial{}, 
			GeomTagToBoundary{}, 
			GeomTagToInteriorConditions{} }, 
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
			GeomTagToMaterial{},
			GeomTagToBoundary{},
			GeomTagToInteriorConditions{}),
		EvolutionOptions{})
	);
}

TEST_F(EigenvalueEstimatorTest, printMatrix_PEC)
{
	EigenvalueEstimator ev(
		fesQuad,
		Model{
			mQuad,
			GeomTagToMaterial{},
			GeomTagToBoundary{},
			GeomTagToInteriorConditions{} },
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
			GeomTagToMaterial{},
			GeomTagToBoundary{ {1,BdrCond::SMA} },
			GeomTagToInteriorConditions{} },
		EvolutionOptions{});

	auto mat{ ev.getElementMatrix() };

	ASSERT_GE(mat.rows(), 0);
	ASSERT_GE(mat.cols(), 0);

	std::cout << mat << std::endl;
	std::cout << std::flush;

}

TEST_F(EigenvalueEstimatorTest, comparePECandSMAconditions)
{
	auto evPEC{EigenvalueEstimator(fesQuad, Model{ mQuad, GeomTagToMaterial{}, GeomTagToBoundary{}, GeomTagToInteriorConditions{} }, EvolutionOptions{}) };
	auto att_to_bdr{ GeomTagToBoundary{{1,BdrCond::SMA}} };
	auto evSMA{ EigenvalueEstimator(fesQuad, Model{ mQuad, GeomTagToMaterial{}, att_to_bdr, GeomTagToInteriorConditions{} }, EvolutionOptions{}) };

	auto eigenvalsPEC{ evPEC.getElementMatrix().eigenvalues() };
	auto eigenvalsSMA{ evSMA.getElementMatrix().eigenvalues() };

	std::cout << "PEC" << std::endl;
	std::cout << evPEC.getElementMatrix().eigenvalues() << std::endl;
	std::cout << "SMA" << std::endl;
	std::cout << evSMA.getElementMatrix().eigenvalues() << std::endl;
}