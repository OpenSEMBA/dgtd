#include "components/EigenvalueEstimator.h"

using namespace maxwell;
using namespace mfem;

using Solver = maxwell::Solver;

class EigenvalueEstimatorTest : public ::testing::Test {
protected:
	Mesh smTri { Mesh::MakeCartesian2D(1, 1,    Element::TRIANGLE     , true) };
	Mesh smQuad{ Mesh::MakeCartesian2D(1, 1,    Element::QUADRILATERAL, true) };
	Mesh smHexa{ Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON) };
	ParMesh mTri = ParMesh(MPI_COMM_WORLD, smTri);
	ParMesh mQuad = ParMesh(MPI_COMM_WORLD, smQuad);
	ParMesh mHexa = ParMesh(MPI_COMM_WORLD, smHexa);
	DG_FECollection fec2D{ DG_FECollection(1, 2,BasisType::GaussLegendre) };
	DG_FECollection fec3D{ DG_FECollection(1, 3,BasisType::GaussLegendre) };
	ParFiniteElementSpace fesTri  = ParFiniteElementSpace(&mTri,  &fec2D);
	ParFiniteElementSpace fesQuad = ParFiniteElementSpace(&mQuad, &fec2D);
	ParFiniteElementSpace fesHexa = ParFiniteElementSpace(&mHexa, &fec3D);

};

TEST_F(EigenvalueEstimatorTest, verifyNoThrow_Quad)
{
	auto model {Model(mQuad,GeomTagToMaterialInfo(),GeomTagToBoundaryInfo())};
	auto eo {EvolutionOptions{}};
	EXPECT_NO_THROW(EigenvalueEstimator(
		fesQuad, 
		model,
		eo)
	);
}

TEST_F(EigenvalueEstimatorTest, verifyNoThrow_SubMesh_Tri)
{
	mTri.SetAttribute(0, 300);
	Array<int> atts(1); atts[0] = 300;
	auto smTri(ParSubMesh::CreateFromDomain(mTri, atts));
	smTri.SetAttribute(0, 1);
	ParFiniteElementSpace smFesTri{ ParFiniteElementSpace(&smTri, &fec2D) };
	
	auto model { Model(mTri,GeomTagToMaterialInfo(),GeomTagToBoundaryInfo()) };
	auto eo {EvolutionOptions{}};
	EXPECT_NO_THROW(EigenvalueEstimator(
		smFesTri, 
		model,
		eo)
	);
}

TEST_F(EigenvalueEstimatorTest, printMatrix_PEC)
{
	auto model{ Model(mQuad,GeomTagToMaterialInfo(),GeomTagToBoundaryInfo()) };
	auto eo {EvolutionOptions{}};
	EigenvalueEstimator ev(
		fesQuad,
		model,
		eo);

	auto mat{ ev.getElementMatrix() };
	
	ASSERT_GE(mat.rows(), 0);
	ASSERT_GE(mat.cols(), 0);

	std::cout << mat << std::endl;
	std::cout << std::flush;

}

TEST_F(EigenvalueEstimatorTest, printMatrix_SMA)
{
	auto att_to_bdr {GeomTagToBoundary{}};
	att_to_bdr.insert({1,BdrCond::SMA});
	auto model {Model(mQuad, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo(att_to_bdr, GeomTagToInteriorBoundary()))};
	auto eo {EvolutionOptions{}};
	EigenvalueEstimator ev(
		fesQuad,
		model,
		eo);

	auto mat{ ev.getElementMatrix() };

	ASSERT_GE(mat.rows(), 0);
	ASSERT_GE(mat.cols(), 0);

	std::cout << mat << std::endl;
	std::cout << std::flush;

}

TEST_F(EigenvalueEstimatorTest, comparePECandSMAconditions)
{
	auto eo {EvolutionOptions{}};
	auto model{ Model{mQuad, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo()} };
	auto evPEC{EigenvalueEstimator(fesQuad, model, eo) };
	auto att_to_bdr{ GeomTagToBoundary{}};
	att_to_bdr.insert({1,BdrCond::SMA});
	auto smaModel{ Model{mQuad, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo(att_to_bdr, GeomTagToInteriorBoundary{})} };
	auto evSMA{ EigenvalueEstimator(fesQuad, smaModel, eo) };

	auto eigenvalsPEC{ evPEC.getElementMatrix().eigenvalues() };
	auto eigenvalsSMA{ evSMA.getElementMatrix().eigenvalues() };

	std::cout << "PEC" << std::endl;
	std::cout << evPEC.getElementMatrix().eigenvalues() << std::endl;
	std::cout << "SMA" << std::endl;
	std::cout << evSMA.getElementMatrix().eigenvalues() << std::endl;
}
