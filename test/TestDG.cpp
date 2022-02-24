#include "gtest/gtest.h"
#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "../../mfemlib4.3/miniapps/common/mfem-common.hpp"
#include "../../mfemlib4.3/miniapps/common/fem_extras.hpp"
#include "../../mfemlib4.3/miniapps/common/fem_extras.cpp"
#include <vector>

using namespace mfem;

namespace HelperFunctions {

	Array<int> getH1LexOrder(const H1_FECollection* fec) 
	{
		auto* fe = fec->FiniteElementForGeometry(Geometry::SEGMENT);
		const NodalFiniteElement* nodal_fe =
			dynamic_cast<const NodalFiniteElement*>(fe);
		Array<int> lexOrder = nodal_fe->GetLexicographicOrdering();
		return lexOrder;
	}

	SparseMatrix operatorToSparseMatrix(const Operator* op) 
	{

		int width = op->Width();
		int height = op->Height();
		SparseMatrix res(height, height);
		Vector x(width), y(height);

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

	void compareH1AndDGMassMatrixes(int& order, Mesh& mesh, const int& basis) 
	{

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

	Mesh buildCartesianMeshForOneElement(const int& dimension, const Element::Type& element) 
	{
		switch (dimension) {
		case 1:
			switch (element) {
			case Element::SEGMENT:
				return Mesh::MakeCartesian1D(1);
				break;
			default:
				throw std::exception("1-Dimensional meshes can only be SEGMENT based.");
			}
		case 2:
			switch (element) {
			case Element::TRIANGLE:
				return Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
				break;
			case Element::QUADRILATERAL:
				return Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL);
				break;
			default:
				throw std::exception("2-Dimensional meshes can only be TRIANGLE or QUADRILATERAL based.");
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
				throw std::exception("3-Dimensional meshes can only be HEXAEDRON, WEDGE or TETRAHEDRON based.");
			}
		default:
			throw std::exception("Dimension must be 2 or 3 with Element argument. Or dimension 1, which ignores element.");
		} 
	}

	double linearFunction(const Vector& pos)
	{
		double normalizedPos;
		double leftBoundary = 0.0, rightBoundary = 1.0;
		double length = rightBoundary - leftBoundary;
		normalizedPos = (pos[0] - leftBoundary) / length;

		return 2*normalizedPos;
	}

	void setInitialCondition(GridFunction& sol,std::function<double(const Vector&)> f) 
	{
		sol.ProjectCoefficient(FunctionCoefficient(f));
	}

	void SaveData(GridFunction& gf, const char* filename) {
		gf.Save(filename);
	}

}

class DG : public ::testing::Test {
};
TEST_F(DG, checkDataValueOutsideNodesForOneElementMeshes)
{
	const int dimension = 1;
	const int order = 1;
	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(1,Element::SEGMENT);
	auto fecDG = new DG_FECollection(order, dimension, BasisType::GaussLegendre);
	auto* fesDG = new FiniteElementSpace(&mesh, fecDG);

	GridFunction solution;
	solution.SetSpace(fesDG);
	solution.ProjectCoefficient(FunctionCoefficient(HelperFunctions::linearFunction));
	IntegrationPoint integPoint;
	for (double xVal = 0.0; xVal <= 1; xVal = xVal + 0.1) {
		integPoint.Set(xVal, 0.0, 0.0, 0.0);
		double interpolatedPoint = solution.GetValue(0, integPoint);
		EXPECT_NEAR(xVal * 2, interpolatedPoint,1e-10);
	}
}
TEST_F(DG, checkMassMatrix)
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

TEST_F(DG, checkMassMatrixIsSameForH1andDG)
{
	const int maxOrder = 10;
	int order = 1;
	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(2, Element::QUADRILATERAL);

	for (order; order < maxOrder; order++) {

		ASSERT_EQ(1, mesh.GetNE());

		HelperFunctions::compareH1AndDGMassMatrixes(order, mesh, BasisType::ClosedUniform);
	}
}
TEST_F(DG, checkStiffnessMatrix)
{
	int order = 2;
	const int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new DG_FECollection(order, dimension,BasisType::GaussLobatto);
	fes = new FiniteElementSpace(&mesh, fec);
	//Coefficient* one = new ConstantCoefficient(1.0);

	BilinearForm stiffnessMatrix(fes);
	stiffnessMatrix.AddDomainIntegrator(
			new DerivativeIntegrator(
				*(new ConstantCoefficient(1.0)),0));
	stiffnessMatrix.Assemble();
	stiffnessMatrix.Finalize();

	auto stiffnessSparse = stiffnessMatrix.SpMat();
	auto stiffnessDense = stiffnessSparse.ToDenseMatrix();

	stiffnessDense->Print(std::cout);
	std::cout << std::endl;

	switch (order) {
	case 1:
		EXPECT_NEAR(-0.5, stiffnessDense->Elem(0, 0), 1e-3);
		EXPECT_NEAR(0.5, stiffnessDense->Elem(0, 1), 1e-3);
		EXPECT_NEAR(-0.5, stiffnessDense->Elem(1, 0), 1e-3);
		EXPECT_NEAR(0.5, stiffnessDense->Elem(1, 1), 1e-3);
		break;
	case 2:
		EXPECT_NEAR(-5e-1, stiffnessDense->Elem(0, 0), 1e-3);
		EXPECT_NEAR(6.6667e-1 ,stiffnessDense->Elem(0, 1), 1e-3);
		EXPECT_NEAR(-1.6667e-1 ,stiffnessDense->Elem(0, 2), 1e-3);
		EXPECT_NEAR(-6.6667e-1,stiffnessDense->Elem(1, 0), 1e-3);
		EXPECT_NEAR(0.0 ,stiffnessDense->Elem(1, 1), 1e-3);
		EXPECT_NEAR(6.6667e-1,stiffnessDense->Elem(1, 2), 1e-3);
		EXPECT_NEAR(1.6667e-1,stiffnessDense->Elem(2, 0), 1e-3);
		EXPECT_NEAR(-6.6667e-1,stiffnessDense->Elem(2, 1), 1e-3);
		EXPECT_NEAR(5e-1 ,stiffnessDense->Elem(2, 2), 1e-3);
		break;
	}

}

TEST_F(DG, checkFluxOperators)
{
	int order = 1;
	const int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new DG_FECollection(order, dimension, BasisType::GaussLobatto);
	fes = new FiniteElementSpace(&mesh, fec);

	BilinearForm fluxForm(fes);
}
TEST_F(DG, visualizeGLVISDataForBasisFunctionNodes)
{
	const int dimension = 1;
	const int order = 1;

	char vishost[] = "localhost";
	int  visport = 19916;

	struct VisWinLayout
	{
		int nx;
		int ny;
		int w;
		int h;
	};

	VisWinLayout vwl;
	vwl.nx = 5;
	vwl.ny = 3;
	vwl.w = 250;
	vwl.h = 250;


	bool visualization = true;
	int onlySome = -1;

	std::vector<socketstream*> socket;


	Vector nodalVector(order + 1);
	Vector dofVector(order + 1);
	IntegrationPoint integPoint;
	Array<int> vdofs;

	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(1, Element::SEGMENT);
	auto fecDG = new DG_FECollection(order, dimension);
	auto* fesDG = new FiniteElementSpace(&mesh, fecDG);

	int ndof = fesDG->GetVSize();
	fesDG->GetElementVDofs(0, vdofs);

	int offx = vwl.w + 10, offy = vwl.h + 45; // window offsets

	for (unsigned int i = 0; i < socket.size(); i++)
	{
		*socket[i] << "keys q";
		delete socket[i];
	}

	socket.resize(ndof);
	for (int i = 0; i < ndof; i++)
	{
		socket[i] = new socketstream; socket[i]->precision(8);
	}
	GridFunction** solution = new GridFunction * [ndof];

	for (int i = 0; i < ndof; i++) {
		solution[i] = new GridFunction(fesDG);
		*solution[i] = 0.0;
		(*solution[i])(vdofs[i]) = 1.0;
	}

	int stopAt = ndof;
	bool vec = false;

	for (int i = 0; i < stopAt; i++)
	{
		if (i == 0 && onlySome > 0 && onlySome < ndof)
		{
			i = onlySome - 1;
			stopAt = std::min(ndof, onlySome + 9);
		}

		std::ostringstream oss;
		oss << "DoF " << i + 1;
		mfem::common::VisualizeField(*socket[i], vishost, visport, *solution[i], oss.str().c_str(),
			(i % vwl.nx) * offx, ((i / vwl.nx) % vwl.ny) * offy, vwl.w, vwl.h, "aaAc", vec);
	}
}

TEST_F(DG, printGLVISDataForBasisFunctionNodes)
{
	const int dimension = 1;
	const int order = 1;

	Vector nodalVector(order + 1);
	Vector dofVector(order + 1);
	IntegrationPoint integPoint;
	Array<int> vdofs;

	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(1, Element::SEGMENT);
	auto fecDG = new DG_FECollection(order, dimension);
	auto* fesDG = new FiniteElementSpace(&mesh, fecDG);

	int ndof = fesDG->GetVSize();
	fesDG->GetElementVDofs(0, vdofs);

	GridFunction** solution = new GridFunction * [ndof];

	for (int i = 0; i < ndof; i++) {
		solution[i] = new GridFunction(fesDG);
		*solution[i] = 0.0;
		(*solution[i])(vdofs[i]) = 1.0;
		std::string stringName = "L2_O" + std::to_string(order) + "_SEG_N" + std::to_string(i) + ".gf";
		const char* filename = stringName.c_str();
		HelperFunctions::SaveData(*solution[i], filename);
	}
	HelperFunctions::SaveData(**solution, "save.gf");
	mesh.Save("mesh.mesh");
}