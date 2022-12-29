#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>

using namespace mfem;

class TestFiniteElementSpace : public ::testing::Test {
protected:

	typedef std::size_t Direction;

	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian1D(1);
		fec_ = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES1D(
		const int order,
		const int elements = 1)
	{
		mesh_ = Mesh::MakeCartesian1D(elements);
		fec_ = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	}

	void setFES2D(
		const int order,
		const int xElem = 1,
		const int yElem = 1)
	{
		mesh_ = Mesh::MakeCartesian2D(xElem, yElem, Element::Type::QUADRILATERAL);
		fec_ = std::make_unique<DG_FECollection>(order, mesh_.Dimension(), BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;
	
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

	std::unique_ptr<SparseMatrix> rotateMatrixLexico(BilinearForm& matrix)
	{
		const Operator* rotatorOperator = matrix.FESpace()->GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC);
		const SparseMatrix rotatorMatrix = operatorToSparseMatrix(rotatorOperator);
		const SparseMatrix matrixSparse = matrix.SpMat();
		SparseMatrix* res;
		{
			auto aux = Mult(matrixSparse, rotatorMatrix);
			res = TransposeMult(rotatorMatrix, *aux);
		}
		return std::unique_ptr<SparseMatrix>(res);
	}

	static std::string getFilename(const std::string fn)
	{
		return "./testData/" + fn;
	}

	void SaveData(GridFunction& gf, const char* filename) {
		gf.Save(filename);
	}
};

TEST_F(TestFiniteElementSpace, MassMatrixIsSameForH1andDG)
{
	for (int order = 1; order < 5; order++) {

		setFES2D(order);

		std::unique_ptr<FiniteElementCollection> fecH1 = std::make_unique<H1_FECollection>(order, mesh_.Dimension(), BasisType::GaussLobatto);
		std::unique_ptr<FiniteElementSpace> fesH1 = std::make_unique<FiniteElementSpace>(&mesh_, fecH1.get());
		BilinearForm massMatrixH1(fesH1.get());
		massMatrixH1.AddDomainIntegrator(new MassIntegrator);
		massMatrixH1.Assemble();
		massMatrixH1.Finalize();

		BilinearForm massMatrixDG(fes_.get());
		massMatrixDG.AddDomainIntegrator(new MassIntegrator);
		massMatrixDG.Assemble();
		massMatrixDG.Finalize();

		for (int i = 0; i < massMatrixDG.SpMat().NumRows(); i++) {
			for (int j = 0; j < massMatrixDG.SpMat().NumCols(); j++) {
				EXPECT_NEAR(rotateMatrixLexico(massMatrixH1)->Elem(i, j), massMatrixDG.SpMat().Elem(i, j), 1e-5);
			}
		}
	}
}
TEST_F(TestFiniteElementSpace, KOperators)
{
	/* The objetive of this test is to check the construction of the bilinear form 
	  for a single element in 1D.

	  Firstly, we build a matrix which represent the total bilinear form of the problem (K).
	  Secondly, we build the stifness (S) and flux (F) matrix and convert all the matrix to 
	  a dense for comparing purporse.
	  
	  Finally, we compare the elements of the initial bilinear form (K) and the sum of the 
	  elements of the the stiffness (S) and flux (F) matrix. */ 

	setFES1D(2);

	ConstantCoefficient one{ 1.0 };
	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	BilinearForm kMat(fes_.get());
	kMat.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, 0)));
	kMat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	kMat.AddBdrFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	kMat.Assemble();
	kMat.Finalize();

	BilinearForm sMat(fes_.get());
	sMat.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, 0)));
	sMat.Assemble();
	sMat.Finalize();

	BilinearForm fMat(fes_.get());
	fMat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	fMat.AddBdrFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	fMat.Assemble();
	fMat.Finalize();

	for (int i = 0; i < kMat.SpMat().ToDenseMatrix()->NumRows(); i++) {
		for (int j = 0; j < kMat.SpMat().ToDenseMatrix()->NumCols(); j++) {
			EXPECT_NEAR(kMat.SpMat().ToDenseMatrix()->Elem(i, j), sMat.SpMat().ToDenseMatrix()->Elem(i, j) + fMat.SpMat().ToDenseMatrix()->Elem(i, j), 1e-3);
		}
	}
}
TEST_F(TestFiniteElementSpace, MeshBoundaries)
{
	/*In this test we aim to compare the boundary DoFs for a mesh generated through
	the Mesh class constructors for a 2DCartesian and 'square3x3.mesh', a handcrafted mesh.

	First, we declare and initialise the variables we will require in order to create a
	FiniteElementSpace object, that is a FiniteElementCollection and a Mesh. In this we create
	two FES, one for the Cartesian2D mesh (Auto) and one for the handcrafted mesh (Manual).

	We then extract the Boundary DoFs from each of the FES.

	Lastly, we compare both of the Arrays where the DoFs are stored, to confirm that when it comes
	to boundaries of this nature, the Cartesian2D constructor works equally to the handcrafted mesh.
	*/

	int order = 1;
	int dimension = 2;
	int nx = 3; 
	int ny = 3; 
	bool generateEdges = true;
	Mesh meshAuto{ Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, generateEdges) };
	Mesh meshManual(*getFilename("square3x3.mesh").c_str(), 1, 1);

	auto fec = new DG_FECollection(order, dimension, BasisType::GaussLegendre);
	auto fesAuto = new FiniteElementSpace(&meshAuto, fec);
	auto fesManual = new FiniteElementSpace(&meshManual, fec);

	Array<int> ess_tdof_list_auto;
	if (meshAuto.bdr_attributes.Size())
	{
		Array<int> ess_bdr_auto(meshAuto.bdr_attributes.Max());
		ess_bdr_auto = 1;
		fesAuto->GetEssentialTrueDofs(ess_bdr_auto, ess_tdof_list_auto);
	}

	Array<int> ess_tdof_list_manual;
	if (meshManual.bdr_attributes.Size())
	{
		Array<int> ess_bdr_manual(meshManual.bdr_attributes.Max());
		ess_bdr_manual = 1;
		fesManual->GetEssentialTrueDofs(ess_bdr_manual, ess_tdof_list_manual);
	}

	EXPECT_EQ(ess_tdof_list_auto, ess_tdof_list_manual);
}
TEST_F(TestFiniteElementSpace, printGLVISDataForBasisFunctionNodes)
{
	/*This test creates files for the Basis Functions, for later visualization 
	through GLVIS.
	
	First, the basic variables and objects to create a FiniteElementSpace are
	declared.

	Then, the number of VDoFs are extracted from the fes.

	Lastly, for each DoF, a 1.0 is assigned to only one of the nodes, while the rest
	are left at 0.0, then a file is written through the GridFunction SaveData function
	for each of the Basis Functions. The mesh too, is saved as a file.
	*/
	
	const int dimension = 1;
	const int order = 2;

	Vector nodalVector(order + 1);
	Vector dofVector(order + 1);
	Array<int> vdofs;

	Mesh mesh = Mesh::MakeCartesian1D(1);
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
		SaveData(*solution[i], stringName.c_str());
	}
	SaveData(**solution, "save.gf");
	mesh.Save("mesh.mesh");
}


