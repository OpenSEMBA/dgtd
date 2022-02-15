#include "gtest/gtest.h"

#include "mfem.hpp"

using namespace mfem;

namespace HelperFunctions {

	mfem::Array<int> getH1LexOrder(const mfem::H1_FECollection* fec) {
		auto* fe = fec->FiniteElementForGeometry(mfem::Geometry::SEGMENT);
		const mfem::NodalFiniteElement* nodal_fe =
			dynamic_cast<const mfem::NodalFiniteElement*>(fe);
		mfem::Array<int> lexOrder = nodal_fe->GetLexicographicOrdering();
		return lexOrder;
	}

	mfem::SparseMatrix operatorToSparseMatrix(const mfem::Operator* op) {

		int width = op->Width();
		int height = op->Height();
		mfem::SparseMatrix res(height, height);
		mfem::Vector x(width), y(height);

		x = 0.0;

		for (int i = 0; i < width; i++)
		{
			x(i) = 1.0;
			op->Mult(x, y);
			for (int j = 0; j < height; j++)
			{
				if (y(j) != 0.0)
				{
					res.Add(i, j, y[j]);
				}
			}
			x(i) = 0.0;
		}

		res.Finalize();
		return res;
	}
	void compareH1AndDGMassMatrixes(int& order, Mesh& mesh, const int& basis) {

		std::cout << "Checking order: " << order << std::endl;

		auto fecH1 = new H1_FECollection(order, mesh.Dimension(), basis);
		auto fecDG = new DG_FECollection(order, mesh.Dimension(), basis);

		FiniteElementSpace* fesH1 = new FiniteElementSpace(
			&mesh, fecH1);
		FiniteElementSpace* fesDG = new FiniteElementSpace(
			&mesh, fecDG);

		const Operator* rotatorH1 = fesH1->GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC);

		const SparseMatrix rotatorMatrix = HelperFunctions::operatorToSparseMatrix(rotatorH1);
		rotatorMatrix.Finalized();

		rotatorMatrix.Print(std::cout);
		std::cout << std::endl;

		BilinearForm massMatrixH1(fesH1);
		massMatrixH1.AddDomainIntegrator(new MassIntegrator);
		massMatrixH1.Assemble();
		massMatrixH1.Finalize();
		BilinearForm massMatrixDG(fesDG);
		massMatrixDG.AddDomainIntegrator(new MassIntegrator);
		massMatrixDG.Assemble();
		massMatrixDG.Finalize();

		const SparseMatrix massMatrixH1Sparse = massMatrixH1.SpMat();
		massMatrixH1Sparse.Finalized();		

		SparseMatrix* rotatedMassMatrixH1Sparse;
		{
			auto aux = Mult(massMatrixH1Sparse, rotatorMatrix);
			rotatedMassMatrixH1Sparse = TransposeMult(rotatorMatrix, *aux);			
		}
		rotatedMassMatrixH1Sparse->Finalized();

		const SparseMatrix massMatrixDGSparse = massMatrixDG.SpMat();
		massMatrixDGSparse.Finalized();
		
		ASSERT_EQ(rotatedMassMatrixH1Sparse->NumRows(), massMatrixDGSparse.NumRows());
		ASSERT_EQ(rotatedMassMatrixH1Sparse->NumCols(), massMatrixDGSparse.NumCols());

		for (std::size_t i = 0; i < massMatrixDGSparse.NumRows(); i++) {
			for (std::size_t j = 0; j < massMatrixDGSparse.NumCols(); j++) {
				EXPECT_NEAR(rotatedMassMatrixH1Sparse->Elem(i, j), massMatrixDGSparse.Elem(i, j), 1e-5);
			}
		}
	}

	Mesh buildCartesianMeshForOneElement(const int& dimension, const Element::Type& element) {
		switch (dimension) {
		case 1:
			return Mesh::MakeCartesian1D(1);
			break;
		case 2:
			switch (element) {
			case Element::TRIANGLE:
				return Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
				break;
			case Element::QUADRILATERAL:
				return Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL);
				break;
			default:
				throw std::exception("Incorrect combination of dimension and/or element type.");
			}
		case 3:
			switch (element) {
			case Element::HEXAHEDRON:
				return Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON);
				break;
			case Element::WEDGE:
				return Mesh::MakeCartesian3D(1, 1, 1, Element::WEDGE);
				break;
			case Element::TETRAHEDRON:
				return Mesh::MakeCartesian3D(1, 1, 1, Element::TETRAHEDRON);
				break;
			default:
				throw std::exception("Incorrect combination of dimension and/or element type.");
			}
		default:
			throw std::exception("Dimension must be 2 or 3 with Element argument. Or dimension 1, which ignores element.");
		} 
	}
}

TEST(DG, checkMassMatrix)
{
	int order = 1;
	int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new H1_FECollection(order, dimension);
	fes = new FiniteElementSpace(&mesh, fec);

	BilinearForm massMatrix(fes);
	massMatrix.AddDomainIntegrator(new MassIntegrator);
	massMatrix.Assemble();
	massMatrix.Finalize();

	EXPECT_NEAR(2.0 / 6.0, massMatrix(0, 0), 1e-3);
	EXPECT_NEAR(1.0 / 6.0, massMatrix(0, 1), 1e-3);
	EXPECT_NEAR(1.0 / 6.0, massMatrix(1, 0), 1e-3);
	EXPECT_NEAR(2.0 / 6.0, massMatrix(1, 1), 1e-3);

}

TEST(DG, checkMassMatrixIsSameForH1andDG)
{
	const int maxOrder = 10;
	int order = 1;
	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(2, Element::QUADRILATERAL);

	for (order; order < maxOrder; order++) {

		ASSERT_EQ(1, mesh.GetNE());

		HelperFunctions::compareH1AndDGMassMatrixes(order, mesh, BasisType::ClosedUniform);
	}
}

TEST(DG, checkDataValueOutsideNodes)
{
	const int dimension = 1;
	const int order = 1;
	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(2,Element::QUADRILATERAL);

	auto fecH1 = new H1_FECollection(order, dimension, BasisType::GaussLegendre);
	auto* fesH1 = new FiniteElementSpace(&mesh, fecH1);

	BilinearForm massMatrixH1(fesH1);
	massMatrixH1.AddDomainIntegrator(new MassIntegrator);
	massMatrixH1.Assemble();
	massMatrixH1.Finalize();





}