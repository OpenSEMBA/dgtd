#include "gtest/gtest.h"
#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "../src/maxwell/BilinearIntegrators.h"
#include "../src/maxwell/Types.h"

#include <vector>
#include <Eigen/Dense>
#include <maxwell/Model.h>

using namespace mfem;

struct FluxCoefficient {
	double alpha;
	double beta;
};

namespace HelperFunctions {

	Mesh makeTwoAttributeCartesianMesh1D(
		const int& refTimes = 0)
	{
		Mesh res = Mesh::MakeCartesian1D(2);
		res.SetAttribute(0, 1);
		res.SetAttribute(1, 2);

		for (int i = 0; i < refTimes; i++) {
			res.UniformRefinement();
		}

		return res;
	}

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

	std::unique_ptr<SparseMatrix> rotateMatrixLexico(BilinearForm& matrix)
	{
		const Operator* rotatorOperator = matrix.FESpace()->GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC);
		const SparseMatrix rotatorMatrix = HelperFunctions::operatorToSparseMatrix(rotatorOperator);
		const SparseMatrix matrixSparse = matrix.SpMat();
		SparseMatrix* res;
		{
			auto aux = Mult(matrixSparse, rotatorMatrix);
			res = TransposeMult(rotatorMatrix, *aux);
		}
		return std::unique_ptr<SparseMatrix>(res);
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

	void checkDenseMatrixSubtractIsValueForAllElem(
		const double val, 
		std::unique_ptr<DenseMatrix> m1, 
		std::unique_ptr<DenseMatrix> m2)
	{
		for (int i = 0; i < m1->Width(); i++) {
			for (int j = 0; j < m1->Height(); j++) {
				EXPECT_NEAR(0.0, m1->Elem(i,j) - m2->Elem(i,j), 1e-3);
			}
		}
	}

	Eigen::Matrix<double, 27, 27> build3DOneElementDMatrix()
	{
		auto res = Eigen::Matrix<double, 27, 27>();
		res.setZero();
		auto blockMat = Eigen::Matrix3d{
				{-1.5, 2.0,-0.5},
				{-0.5, 0.0, 0.5},
				{ 0.5,-2.0, 1.5} };
		for (int i = 0; i < res.cols(); i += 3) {
			res.block<3, 3>(i,i) = blockMat;
		}
		return res;
	}

	
}

class TestMFEMFunctionality : public ::testing::Test {
protected:
	typedef std::size_t Direction;
	
};

TEST_F(TestMFEMFunctionality, checkMassMatrixIsSameForH1andDG)
{
	/*This test compares the mass matrices for H1 and DG spaces for a mesh with a single element

	First, a mesh is built with buildCartesianMeshForOneElement() which ensures a mesh
	with a single element will be used. Then, for different orders, FiniteElementCollection and
	FiniteElementSpace will be created for both H1 and DG, along with a BilinearForm containing
	the Mass Matrix for each one of them.

	Then, the BilinearForm for the Mass Matrix for H1 will have its inner Sparse Matrix rotated
	by applying a lexicographic rotation operator, as spaces for H1 and DG are ordered differently,
	and be returned as a Sparse Matrix. The DG Bilinear Form will just have its Sparse Matrix
	extracted.

	Lastly, assertions will be made to ensure the number of rows and columns are the same
	for both matrices. To finalise, each element will be compared for the same position in the
	matrices.*/

	const int maxOrder = 5;
	int order = 1;
	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(2, Element::QUADRILATERAL);

	for (order; order < maxOrder; order++) {

		ASSERT_EQ(1, mesh.GetNE());

		std::cout << "Checking order: " << order << std::endl;

		auto fecH1 = new H1_FECollection(order, mesh.Dimension(), BasisType::ClosedUniform);
		FiniteElementSpace* fesH1 = new FiniteElementSpace(&mesh, fecH1);
		BilinearForm massMatrixH1(fesH1);
		massMatrixH1.AddDomainIntegrator(new MassIntegrator);
		massMatrixH1.Assemble();
		massMatrixH1.Finalize();

		auto fecDG = new DG_FECollection(order, mesh.Dimension(), BasisType::ClosedUniform);
		FiniteElementSpace* fesDG = new FiniteElementSpace(&mesh, fecDG);
		BilinearForm massMatrixDG(fesDG);
		massMatrixDG.AddDomainIntegrator(new MassIntegrator);
		massMatrixDG.Assemble();
		massMatrixDG.Finalize();

		auto rotatedMassMatrixH1Sparse = HelperFunctions::rotateMatrixLexico(massMatrixH1);
		auto massMatrixDGSparse = massMatrixDG.SpMat();

		ASSERT_EQ(rotatedMassMatrixH1Sparse->NumRows(), massMatrixDGSparse.NumRows());
		ASSERT_EQ(rotatedMassMatrixH1Sparse->NumCols(), massMatrixDGSparse.NumCols());

		for (int i = 0; i < massMatrixDGSparse.NumRows(); i++) {
			for (int j = 0; j < massMatrixDGSparse.NumCols(); j++) {
				EXPECT_NEAR(rotatedMassMatrixH1Sparse->Elem(i, j), massMatrixDGSparse.Elem(i, j), 1e-5);
			}
		}
	}
}

//TEST_F(Auxiliary, checkDOperator3D)
//{
//
//	Mesh mesh = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::QUADRILATERAL);
//	FiniteElementCollection* fec = new DG_FECollection(2, 3, BasisType::GaussLobatto);
//	FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);
//
//	auto mInvEigen = HelperFunctions::buildInverseMassMatrixEigen(fes);
//	auto oldEigen = HelperFunctions::buildStiffnessMatrixEigen(fes);
//	auto MISCalcOld = 0.5 * mInvEigen * oldEigen;
//	auto expMat = HelperFunctions::build3DOneElementDMatrix();
//
//	std::cout << expMat << std::endl;
//
//	EXPECT_TRUE(MISCalcOld.isApprox(expMat, 1e-1));
//}
TEST_F(TestMFEMFunctionality, checkTwoAttributeMesh)
{
	/*The purpose of this test is to check the makeTwoAttributeCartesianMesh1D(const int& refTimes)
	function.

	First, an integer is declared for the number of times we wish to refine the mesh, then a mesh is
	constructed with two elements, left and right hand sides, setting the following attributes.

	|------LHS------|------RHS------|

	|##ATTRIBUTE 1##|##ATTRIBUTE 2##|

	Once the mesh is refined, it is returned, then we compare if the expected number of elements is
	true for the actual elements in the mesh.

	Then, we consider how the mesh will perform its uniform refinement, and we declare that the
	LHS elements with Attribute one will be Even index elements (starting at 0), and the RHS
	elements with Attribute 2 will be Uneven index elements (starting at 1).*/

	const int refTimes = 3;
	Mesh mesh = HelperFunctions::makeTwoAttributeCartesianMesh1D(refTimes);

	EXPECT_EQ(pow(2, refTimes + 1), mesh.GetNE());
	for (int i = 0; i < mesh.GetNE(); i++) {
		if (i % 2 == 0) {
			EXPECT_EQ(1, mesh.GetAttribute(i));
		}
		else {
			EXPECT_EQ(2, mesh.GetAttribute(i));
		}
	}
}
TEST_F(TestMFEMFunctionality, checkKOperators)
{
	/* The objetive of this test is to check the construction of the bilinear form 
	  for a single element in 1D.

	  Firstly, we build a matrix which represent the total bilinear form of the problem (K).
	  Secondly, we build the stifness (S) and flux (F) matrix and convert all the matrix to 
	  a dense for comparing purporse.
	  
	  Finally, we compare the elements of the initial bilinear form (K) and the sum of the 
	  elements of the the stiffness (S) and flux (F) matrix. */ 

	Mesh mesh = Mesh::MakeCartesian1D(1);
	FiniteElementCollection* fec = new DG_FECollection(2, 1, BasisType::GaussLobatto);
	FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);

	ConstantCoefficient one{ 1.0 };
	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	BilinearForm kMat(fes);
	kMat.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, 0)));
	kMat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	kMat.AddBdrFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	kMat.Assemble();
	kMat.Finalize();

	BilinearForm sMat(fes);
	sMat.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, 0)));
	sMat.Assemble();
	sMat.Finalize();

	BilinearForm fMat(fes);
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

TEST_F(TestMFEMFunctionality, printGLVISDataForBasisFunctionNodes)
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
TEST_F(TestMFEMFunctionality, checkDataValueOutsideNodesForOneElementMeshes)
{
	/* The purpose of this test is to ensure we can extract data from a GridFunction,
	even if the point we're trying to obtain it at is not necessarily a DoF or node.
	
	First, the basic process to declare and initialise a FiniteElementSpace is done,
	this means variables such as dimension, order, creating a mesh, a FEC and a finally,
	the FES.
	
	A GridFunction is then created and assigned the FES. A function is projected in the
	GridFunction, which is a linear function with a slope of 2.
	
	Lastly, an IntegrationPoint is constructed, which we will use to obtain the values
	from the GridFunction at any point we want. As the slope of the line is 2, we expect
	the values to be 2 times the xVal.*/
	
	const int dimension = 1;
	const int order = 1;
	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(1,Element::SEGMENT);
	auto fecDG = new DG_FECollection(order, dimension, BasisType::GaussLegendre);
	auto* fesDG = new FiniteElementSpace(&mesh, fecDG);

	GridFunction solution(fesDG);
	solution.ProjectCoefficient(FunctionCoefficient(HelperFunctions::linearFunction));
	IntegrationPoint integPoint;
	for (double xVal = 0.0; xVal <= 1; xVal = xVal + 0.1) {
		integPoint.Set(xVal, 0.0, 0.0, 0.0);
		double interpolatedPoint = solution.GetValue(0, integPoint);
		EXPECT_NEAR(xVal * 2, interpolatedPoint,1e-10);
	}
}
TEST_F(TestMFEMFunctionality, findPointsTest)
{
	Mesh mesh = Mesh::MakeCartesian3D(2, 4, 6, Element::Type::HEXAHEDRON, 2.0, 4.0, 6.0);
	DenseMatrix pointMat({ { 0.2,0.4,0.6 },{1.5, 3.5, 5.5},{0.25, 1.25, 3.75},{2.0, 4.0, 6.0} });
	Array<int> elArray;
	Array<IntegrationPoint> ipArray;
	std::vector<double> expVals({ 0.2,0.4,0.6 });

	pointMat.Transpose();
	mesh.FindPoints(pointMat, elArray, ipArray);

	EXPECT_EQ(3, pointMat.Height());
	EXPECT_EQ(4, pointMat.Width());
	EXPECT_EQ(0, elArray[0]);
	EXPECT_EQ(expVals[0], ipArray[0].x);
	EXPECT_EQ(expVals[1], ipArray[0].y);
	EXPECT_EQ(expVals[2], ipArray[0].z);

}
