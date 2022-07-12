#include "gtest/gtest.h"
#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "../src/maxwell/BilinearIntegrators.h"
#include "../src/maxwell/Types.h"

#include <vector>

using namespace mfem;

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

	SparseMatrix* rotateMatrixLexico(BilinearForm& matrix)
	{
		const Operator* rotatorOperator = matrix.FESpace()->GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC);
		const SparseMatrix rotatorMatrix = HelperFunctions::operatorToSparseMatrix(rotatorOperator);
		const SparseMatrix matrixSparse = matrix.SpMat();
		SparseMatrix* res;
		{
			auto aux = Mult(matrixSparse, rotatorMatrix);
			res = TransposeMult(rotatorMatrix, *aux);
		}
		return res;
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

	FiniteElementSpace* buildBilinearFormWith1DCartesianMesh(const int elements, const int order) {
		Mesh mesh = Mesh::MakeCartesian1D(elements);
		FiniteElementCollection* fec = new DG_FECollection(order, mesh.Dimension(), BasisType::GaussLobatto);
		FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);
		return fes;
	}

}

class Auxiliary : public ::testing::Test {
protected:
	typedef std::size_t Direction;
};

TEST_F(Auxiliary, checkDataValueOutsideNodesForOneElementMeshes)
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

TEST_F(Auxiliary, checkVDIM)
{
	Mesh mesh1D = Mesh::MakeCartesian1D(5);
	Mesh mesh2D = Mesh::MakeCartesian2D(5, 5, Element::Type::QUADRILATERAL);
	Mesh mesh3D = Mesh::MakeCartesian3D(5, 5, 5, Element::Type::TETRAHEDRON);

	DG_FECollection fec1D = DG_FECollection(1, mesh1D.Dimension());
	DG_FECollection fec2D = DG_FECollection(1, mesh2D.Dimension());
	DG_FECollection fec3D = DG_FECollection(1, mesh3D.Dimension());

	FiniteElementSpace fes1D = FiniteElementSpace(&mesh1D, &fec1D);
	FiniteElementSpace fes2D = FiniteElementSpace(&mesh2D, &fec2D);
	FiniteElementSpace fes3D = FiniteElementSpace(&mesh3D, &fec3D);

	GridFunction gf1D(&fes1D);
	GridFunction gf2D(&fes2D);
	GridFunction gf3D(&fes3D);

}
TEST_F(Auxiliary, checkMassMatrix)
{
	/*The purpose of this text is to check the values of a Mass Matrix 
	for an order 1, 1D line with a single element.
	
	First, the basic variables and objects to create a FiniteElementSpace are 
	declared.
	
	Then, a BilinearForm is created, in which we add a domain integrator for the
	Mass Matrix, which is MassIntegrator.
	
	Then, we compare the values of the Mass Matrix with those found in Silvester,
	Appendix 3 for a 1D Line Segment.*/
	
	const double tol = 1e-3;
	int order = 1;
	const int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new H1_FECollection(order, dimension);
	fes = new FiniteElementSpace(&mesh, fec);

	BilinearForm massMatrix(fes);
	massMatrix.AddDomainIntegrator(new MassIntegrator);
	massMatrix.Assemble();
	massMatrix.Finalize();

	EXPECT_NEAR(2.0 / 6.0, massMatrix(0, 0), tol);
	EXPECT_NEAR(1.0 / 6.0, massMatrix(0, 1), tol);
	EXPECT_NEAR(1.0 / 6.0, massMatrix(1, 0), tol);
	EXPECT_NEAR(2.0 / 6.0, massMatrix(1, 1), tol);

}

TEST_F(Auxiliary, checkInverseMassMatrix)
{
	/*The purpose of this text is to check the values of a Mass Matrix
	for an order 1, 1D line with a single element.

	First, the basic variables and objects to create a FiniteElementSpace are
	declared.

	Then, a BilinearForm is created, in which we add a domain integrator for the
	Mass Matrix, which is MassIntegrator.

	Then, we compare the values of the Mass Matrix with those found in Silvester,
	Appendix 3 for a 1D Line Segment.*/

	const double tol = 1e-3;
	int order = 1;
	const int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new H1_FECollection(order, dimension);
	fes = new FiniteElementSpace(&mesh, fec);

	BilinearForm massMatrix(fes);
	massMatrix.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	massMatrix.Assemble();
	massMatrix.Finalize();

	EXPECT_NEAR(4.0, massMatrix(0, 0), tol);
	EXPECT_NEAR(-2.0, massMatrix(0, 1), tol);
	EXPECT_NEAR(-2.0, massMatrix(1, 0), tol);
	EXPECT_NEAR(4.0, massMatrix(1, 1), tol);

}

TEST_F(Auxiliary, checkTwoAttributeMesh)
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

TEST_F(Auxiliary, checkMassMatrixIsSameForH1andDG)
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
		FiniteElementSpace* fesH1 = new FiniteElementSpace(
			&mesh, fecH1);
		BilinearForm massMatrixH1(fesH1);
		massMatrixH1.AddDomainIntegrator(new MassIntegrator);
		massMatrixH1.Assemble();
		massMatrixH1.Finalize();

		auto fecDG = new DG_FECollection(order, mesh.Dimension(), BasisType::ClosedUniform);
		FiniteElementSpace* fesDG = new FiniteElementSpace(
			&mesh, fecDG);
		BilinearForm massMatrixDG(fesDG);
		massMatrixDG.AddDomainIntegrator(new MassIntegrator);
		massMatrixDG.Assemble();
		massMatrixDG.Finalize();

		SparseMatrix* rotatedMassMatrixH1Sparse = HelperFunctions::rotateMatrixLexico(massMatrixH1);
		SparseMatrix massMatrixDGSparse = massMatrixDG.SpMat();

		ASSERT_EQ(rotatedMassMatrixH1Sparse->NumRows(), massMatrixDGSparse.NumRows());
		ASSERT_EQ(rotatedMassMatrixH1Sparse->NumCols(), massMatrixDGSparse.NumCols());

		for (int i = 0; i < massMatrixDGSparse.NumRows(); i++) {
			for (int j = 0; j < massMatrixDGSparse.NumCols(); j++) {
				EXPECT_NEAR(rotatedMassMatrixH1Sparse->Elem(i, j), massMatrixDGSparse.Elem(i, j), 1e-5);
			}
		}
	}
}
TEST_F(Auxiliary, checkStiffnessMatrix)
{
	/*The purpose of this text is to check the values of a Stiffness Matrix
	for a 1D line with a single element.

	First, the basic variables and objects to create a FiniteElementSpace are
	declared.

	Then, a BilinearForm is created, in which we add a domain integrator for the
	Stiffness Matrix, which is DerivativeIntegrator. We extract the sparse matrix
	stored inside the BilinearForm, and then convert it to a Dense, for comparing
	purposes.

	Then, we compare the values of the Stiffness Matrix with those calculated
	through Hesthaven's MatLab code for a 1D Line Segment with a single element.*/
	
	const double tol = 1e-3;
	int order = 1;
	const int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new DG_FECollection(order, dimension,BasisType::GaussLobatto);
	fes = new FiniteElementSpace(&mesh, fec);

	BilinearForm stiffnessMatrix(fes);
	stiffnessMatrix.AddDomainIntegrator(
			new DerivativeIntegrator(
				*(new ConstantCoefficient(1.0)),0));
	stiffnessMatrix.Assemble();
	stiffnessMatrix.Finalize();

	auto stiffnessDense = stiffnessMatrix.SpMat().ToDenseMatrix();

	stiffnessDense->Print(std::cout);
	std::cout << std::endl;

	switch (order) {
	case 1:
		EXPECT_NEAR(-0.5, stiffnessDense->Elem(0, 0), tol);
		EXPECT_NEAR(0.5, stiffnessDense->Elem(0, 1), tol);
		EXPECT_NEAR(-0.5, stiffnessDense->Elem(1, 0), tol);
		EXPECT_NEAR(0.5, stiffnessDense->Elem(1, 1), tol);
		break;
	case 2:
		EXPECT_NEAR(-5e-1, stiffnessDense->Elem(0, 0), tol);
		EXPECT_NEAR(6.6667e-1 ,stiffnessDense->Elem(0, 1), tol);
		EXPECT_NEAR(-1.6667e-1 ,stiffnessDense->Elem(0, 2), tol);
		EXPECT_NEAR(-6.6667e-1,stiffnessDense->Elem(1, 0), tol);
		EXPECT_NEAR(0.0 ,stiffnessDense->Elem(1, 1), tol);
		EXPECT_NEAR(6.6667e-1,stiffnessDense->Elem(1, 2), tol);
		EXPECT_NEAR(1.6667e-1,stiffnessDense->Elem(2, 0), tol);
		EXPECT_NEAR(-6.6667e-1,stiffnessDense->Elem(2, 1), tol);
		EXPECT_NEAR(5e-1 ,stiffnessDense->Elem(2, 2), tol);
		break;
	}
}
TEST_F(Auxiliary, checkKOperators)
{
	/* The objetive of this test is to check the construction of the bilinear form 
	  for a single element in 1D.

	  Firstly, we build a matrix which represent the total bilinear form of the problem (K).
	  Secondly, we build the stifness (S) and flux (F) matrix and convert all the matrix to 
	  a dense for comparing purporse.
	  
	  Finally, we compare the elements of the initial bilinear form (K) and the sum of the 
	  elements of the the stiffness (S) and flux (F) matrix. */ 


	int order = 2;
	const int dimension = 1;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	FiniteElementCollection* fec = new DG_FECollection(order, dimension, BasisType::GaussLobatto);
	FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);

	ConstantCoefficient one(1.0);
	double alpha = -1.0;
	double beta = 0.0;
	const Auxiliary::Direction d = 0;
	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	BilinearForm kMat(fes);
	kMat.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, d)));
	kMat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[d], alpha, beta));
	kMat.AddBdrFaceIntegrator(
		new DGTraceIntegrator(n[d], alpha, beta));
	kMat.Assemble();
	kMat.Finalize();

	BilinearForm sMat(fes);
	sMat.AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, d)));
	sMat.Assemble();
	sMat.Finalize();

	BilinearForm fMat(fes);
	fMat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[d], alpha, beta));
	fMat.AddBdrFaceIntegrator(
		new DGTraceIntegrator(n[d], alpha, beta));
	fMat.Assemble();
	fMat.Finalize();

	auto matrixK = kMat.SpMat().ToDenseMatrix();
	auto matrixS = sMat.SpMat().ToDenseMatrix();
	auto matrixF = fMat.SpMat().ToDenseMatrix();

	ASSERT_EQ(matrixK->NumRows(), matrixS->NumRows());
	ASSERT_EQ(matrixK->NumCols(), matrixS->NumCols());
	ASSERT_EQ(matrixK->NumRows(), matrixF->NumRows());
	ASSERT_EQ(matrixK->NumCols(), matrixF->NumCols());

	const double tol = 1e-3;

	for (int i = 0; i < matrixK->NumRows(); i++) {
		for (int j = 0; j < matrixK->NumCols(); j++) {
			EXPECT_NEAR(matrixK->Elem(i, j), matrixS->Elem(i, j) + matrixF->Elem(i, j), tol);
		}
	}
}

TEST_F(Auxiliary, checkDGTraceAverageOnlyMatrix)
{
	/* This test checks the matrix built by the DGTraceIntegrator for
	an integrator that only consists of Average Operators on 1D.
	The order will be changed each step to verify the behaviour
	is consistent with the supposition that only the four central
	elements of the matrix will be different from zero.

	This test also verifies the Lexicographic ordering of the 
	dofs for each element.

	The system will be composed of a different element mesh with a FE
	Collection of DG elements. A FiniteElementSpace is built with
	the previous elements and then an InteriorFaceIntegrator based
	on a DGTraceIntegrator is added, considering a DGTraceIntegrator
	has the form:
	alpha < rho_u (u.n) {v},[w] > + beta < rho_u |u.n| [v],[w] >
	we declare the arguments to be alpha = 1.0, beta = 0.0 in the X
	direction.

	Lastly, a check for the matrix dimensions is performed, then,
	loop checks if the elements in the center of the matrix
	have the correct values for assembling an average {q} = (q- + q+)/2.
	*/
	for (int elements = 2; elements < 5; elements++) {
		for (int order = 2; order < 5; order++) {

			std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
			auto fes = HelperFunctions::buildBilinearFormWith1DCartesianMesh(elements, order);
			auto DGmat = std::make_unique<BilinearForm>(fes);
			DGmat->AddInteriorFaceIntegrator(
				new DGTraceIntegrator(n[0], 1.0, 0.0));
			DGmat->Assemble();
			DGmat->Finalize();

			DenseMatrix* DGDense = DGmat->SpMat().ToDenseMatrix();

			DGDense->PrintMatlab(std::cout);

			EXPECT_EQ((order + 1) * elements, DGDense->Width());
			EXPECT_EQ((order + 1) * elements, DGDense->Height());

			for (int i = 0; i < order; i++) {
				for (int j = 0; j < order; j++) {
					for (int it = 0; it < elements - 2; it++) {
						if (
							i == order + ((order + 1) * it)     && j == order + ((order + 1) * it) ||
							i == order + ((order + 1) * it)     && j == order + ((order + 1) * it) + 1) {
							EXPECT_NEAR(0.5, DGDense->Elem(i, j), 1e-3);
						}
						else if (
							i == order + ((order + 1) * it) + 1 && j == order + ((order + 1) * it) ||
							i == order + ((order + 1) * it) + 1 && j == order + ((order + 1) * it) + 1) {
							EXPECT_NEAR(-0.5, DGDense->Elem(i , j), 1e-3);
						}
						else {
							EXPECT_NEAR(0.0, DGDense->Elem(i , j), 1e-3);
						}
					}
				}
			}
		}
	}
}
TEST_F(Auxiliary, checkDGTraceJumpOnlyMatrix)
{
	/* This test checks the matrix built by the DGTraceIntegrator for
	   an integrator that only consists of Jump Operators on 1D.           
	   The order will be changed each step to verify the behaviour
	   is consistent with the supposition that only the four central
	   elements of the matrix will be different from zero.

		This test also verifies the Lexicographic ordering of the
		dofs for each element.
	   
	   The system will be composed of a different element mesh with a FE
	   Collection of DG elements. A FiniteElementSpace is built with
	   the previous elements and then an InteriorFaceIntegrator based
	   on a DGTraceIntegrator is added, considering a DGTraceIntegrator
	   has the form:
	   alpha < rho_u (u.n) {v},[w] > + beta < rho_u |u.n| [v],[w] >
	   we declare the arguments to be alpha = 0.0, beta = 1.0 in the X
	   direction.
	   
	   Lastly, a check for the matrix dimensions is performed, then, 
	   loop checks if the elements in the center of the matrix
	   have the correct values for assembling a jump [q] = q- - q+.*/

	for (int elements = 2; elements < 5; elements++) {
		for (int order = 1; order < 5; order++) {
		
			std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
			auto fes = HelperFunctions::buildBilinearFormWith1DCartesianMesh(elements, order);
			auto DGmat = std::make_unique<BilinearForm>(fes);
			DGmat->AddInteriorFaceIntegrator(
				new DGTraceIntegrator(n[0], 0.0, 1.0));
			DGmat->Assemble();
			DGmat->Finalize();

			DenseMatrix* DGDense = DGmat->SpMat().ToDenseMatrix();

			EXPECT_EQ((order + 1) * elements, DGDense->Width());
			EXPECT_EQ((order + 1) * elements, DGDense->Height());

			for (int i = 0; i < order; i++) {
				for (int j = 0; j < order; j++) {
					for (int it = 0; it < elements - 2; it++) {
						if (	 
							i == order + ((order + 1) * it)		&& j == order + ((order + 1) * it)	   ||
							i == order + ((order + 1) * it) + 1 && j == order + ((order + 1) * it) + 1) {
							EXPECT_NEAR(1.0, DGDense->Elem(i, j), 1e-3);
						}
						else if (
							i == order + ((order + 1) * it)		&& j == order + ((order + 1) * it) + 1 ||
							i == order + ((order + 1) * it) + 1 && j == order + ((order + 1) * it)) {
							EXPECT_NEAR(-1.0, DGDense->Elem(i, j), 1e-3);
						}
						else {
							EXPECT_NEAR(0.0, DGDense->Elem(i, j), 1e-3);
						}
					}
				}
			}
		}
	}
}

TEST_F(Auxiliary, checkMaxwellDGTraceJumpOnlyMatrix)
{
	/* This test checks the matrix built by the MaxwellDGTraceIntegrator for
	   an integrator that only consists of Jump Operators on 1D.
	   The order will be changed each step to verify the behaviour
	   is consistent with the supposition that only the four central
	   elements of the matrix will be different from zero.

	   This test also verifies the Lexicographic ordering of the
	   dofs for each element.

	   The system will be composed of a different element mesh with a FE
	   Collection of DG elements. A FiniteElementSpace is built with
	   the previous elements and then an InteriorFaceIntegrator based
	   on a MaxwellDGTraceIntegrator is added, considering a 
	   MaxwellDGTraceIntegrator has the form:
	   beta < rho_u n [v],[w] >
	   we declare the arguments to be Dir = X, beta = 1.0 in the X
	   direction.

	   Lastly, a check for the matrix dimensions is performed, then,
	   loop checks if the elements in the center of the matrix
	   have the correct values for assembling a jump [q] = q- - q+.*/

	for (int elements = 2; elements < 5; elements++) {
		for (int order = 1; order < 5; order++) {

			auto fes = HelperFunctions::buildBilinearFormWith1DCartesianMesh(elements, order);
			auto DGmat = std::make_unique<BilinearForm>(fes);
			DGmat->AddInteriorFaceIntegrator(
				new maxwell::MaxwellDGTraceJumpIntegrator(std::vector<maxwell::Direction>{maxwell::Direction::X}, 1.0));
			DGmat->Assemble();
			DGmat->Finalize();

			DenseMatrix* DGDense = DGmat->SpMat().ToDenseMatrix();

			EXPECT_EQ((order + 1) * elements, DGDense->Width());
			EXPECT_EQ((order + 1) * elements, DGDense->Height());

			for (int i = 0; i < order; i++) {
				for (int j = 0; j < order; j++) {
					for (int it = 0; it < elements - 2; it++) {
						if (
							i == order + ((order + 1) * it) && j == order + ((order + 1) * it) ||
							i == order + ((order + 1) * it) + 1 && j == order + ((order + 1) * it) + 1) {
							EXPECT_NEAR(1.0, DGDense->Elem(i, j), 1e-3);
						}
						else if (
							i == order + ((order + 1) * it) && j == order + ((order + 1) * it) + 1 ||
							i == order + ((order + 1) * it) + 1 && j == order + ((order + 1) * it)) {
							EXPECT_NEAR(-1.0, DGDense->Elem(i, j), 1e-3);
						}
						else {
							EXPECT_NEAR(0.0, DGDense->Elem(i, j), 1e-3);
						}
					}
				}
			}
		}
	}
}

//
//TEST_F(Auxiliary, visualizeGLVISDataForBasisFunctionNodes)
//{
//	/*This test aims to show the Basis Functions through GLVIS visualization.
//	
//	This test has to be reformatted to be properly commented.*/
//	
//	const int dimension = 1;
//	const int order = 1;
//
//	char vishost[] = "localhost";
//	int  visport = 19916;
//
//	struct VisWinLayout
//	{
//		int nx;
//		int ny;
//		int w;
//		int h;
//	};
//
//	VisWinLayout vwl;
//	vwl.nx = 5;
//	vwl.ny = 3;
//	vwl.w = 250;
//	vwl.h = 250;
//
//	bool visualization = true;
//	int onlySome = -1;
//
//	std::vector<socketstream*> socket;
//
//	Vector nodalVector(order + 1);
//	Vector dofVector(order + 1);
//	IntegrationPoint integPoint;
//	Array<int> vdofs;
//
//	Mesh mesh = HelperFunctions::buildCartesianMeshForOneElement(1, Element::SEGMENT);
//	auto fecDG = new DG_FECollection(order, dimension);
//	auto* fesDG = new FiniteElementSpace(&mesh, fecDG);
//
//	int ndof = fesDG->GetVSize();
//	fesDG->GetElementVDofs(0, vdofs);
//
//	int offx = vwl.w + 10, offy = vwl.h + 45; // window offsets
//
//	for (unsigned int i = 0; i < socket.size(); i++)
//	{
//		*socket[i] << "keys q";
//		delete socket[i];
//	}
//
//	socket.resize(ndof);
//	for (int i = 0; i < ndof; i++)
//	{
//		socket[i] = new socketstream; socket[i]->precision(8);
//	}
//	GridFunction** solution = new GridFunction * [ndof];
//
//	for (int i = 0; i < ndof; i++) {
//		solution[i] = new GridFunction(fesDG);
//		*solution[i] = 0.0;
//		(*solution[i])(vdofs[i]) = 1.0;
//	}
//
//	int stopAt = ndof;
//	bool vec = false;
//
//	for (int i = 0; i < stopAt; i++)
//	{
//		if (i == 0 && onlySome > 0 && onlySome < ndof)
//		{
//			i = onlySome - 1;
//			stopAt = std::min(ndof, onlySome + 9);
//		}
//
//		std::ostringstream oss;
//		oss << "DoF " << i + 1;
//		mfem::common::VisualizeField(*socket[i], vishost, visport, *solution[i], oss.str().c_str(),
//			(i % vwl.nx) * offx, ((i / vwl.nx) % vwl.ny) * offy, vwl.w, vwl.h, "aaAc", vec);
//	}
//}

TEST_F(Auxiliary, printGLVISDataForBasisFunctionNodes)
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

TEST_F(Auxiliary, findPointsTest)
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
