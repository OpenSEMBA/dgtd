#include "gtest/gtest.h"
#include <math.h>

#include "maxwell/Solver.h"

using namespace maxwell;

namespace AnalyticalFunctions1D {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;

	const double PI = atan(1.0) * 4;

	double gaussianFunction(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 2.0 * (pos[0] - center) /
			            ((meshBoundingBoxMax[0] - meshBoundingBoxMin[0]));
		
		return exp(-20. * pow(normalizedPos, 2));
	}

	double gaussianFunctionHalfWidth(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 4.0 * (pos[0] - center/2) /
			((meshBoundingBoxMax[0] - meshBoundingBoxMin[0]));

		return exp(-20. * pow(normalizedPos, 2));
	}
}

namespace HelperFunctions1D {

	Mesh makeTwoAttributeCartesianMesh1D(const int& refTimes = 0)
	{
		Mesh res = Mesh::MakeCartesian1D(2);
		res.SetAttribute(0, 1);
		res.SetAttribute(1, 2);

		for (int i = 0; i < refTimes; i++) {
			res.UniformRefinement();
		}

		return res;
	}

	void setAttributeIntervalMesh(const int& attVal, const Vector& elIndexes, Mesh& mesh)
	{
		assert(elIndexes[0] < elIndexes[1], "Lower Index bigger than Higher Index.");
		assert(elIndexes[1] <= mesh.GetNE(), "Declared element index bigger than Mesh Number of Elements.");
			for (int i = elIndexes[0] - 1; i < elIndexes[1]; i++) {
				mesh.SetAttribute(i, attVal);
			}
	}

	std::vector<int> mapQuadElementTopLeftVertex(const mfem::Mesh& mesh)
	{
		std::vector<int> res;
		for (int i = 0; i < mesh.GetNE(); i++) {
			mfem::Array<int> meshArrayElement;
			mesh.GetElementVertices(i, meshArrayElement);
			res.push_back(meshArrayElement[0]);
		}

		return res;
	}

	std::vector<std::pair<attribute, Material>> buildAttToMatVec(const std::vector<attribute>& attVec, const std::vector<Material>& matVec)
	{
		std::vector<std::pair<attribute, Material>> res;
		for (int i = 0; i < attVec.size(); i++) {
			res.push_back(std::make_pair(attVec[i], matVec[i]));
		}
		return res;
	}
}
using namespace AnalyticalFunctions1D;

class TestMaxwellSolver1D : public ::testing::Test {
protected:

	Mesh mesh1D = Mesh::MakeCartesian1D(3,1.0);
	Mesh mesh2D = Mesh::MakeCartesian2D(2, 3, Element::Type::QUADRILATERAL, 2.0, 3.0);
	Mesh mesh3D = Mesh::MakeCartesian3D(2, 4, 6, Element::Type::HEXAHEDRON, 2.0, 4.0, 6.0);
	
	Material mat11 = Material(1.0, 1.0); Material mat12 = Material(1.0, 2.0);
	Material mat21 = Material(2.0, 1.0); Material mat22 = Material(2.0, 2.0);
	
	std::vector<attribute> attArrSingle = std::vector<attribute>({ 1 });
	std::vector<attribute> attArrMultiple = std::vector<attribute>({ 1, 2, 3, 4 });
	std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
	std::vector<Material> matArrMultiple = std::vector<Material>({ mat11,mat12,mat21,mat22 });

	std::vector<std::pair<attribute, Material>> attToMatVec = HelperFunctions1D::buildAttToMatVec(attArrSingle, matArrSimple);

	Model testModel = Model(mesh1D, attToMatVec);

	double spread = 2.0;
	double delay = 0.0;
	Direction d = X;
	FieldType ft = E;

	Source testSource = Source(testModel, spread, delay, d, ft);

	Probes defaultProbes;

	Options defaultOptions;

};

TEST_F(TestMaxwellSolver1D, checkTwoAttributeMesh)
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
	Mesh mesh = HelperFunctions1D::makeTwoAttributeCartesianMesh1D(refTimes);

	EXPECT_EQ(pow(2,refTimes + 1), mesh.GetNE());
	for (int i = 0; i < mesh.GetNE(); i++) {
		if (i % 2 == 0) {
			EXPECT_EQ(1, mesh.GetAttribute(i));
		}
		else {
			EXPECT_EQ(2, mesh.GetAttribute(i));
		}
	}
}
TEST_F(TestMaxwellSolver1D, oneDimensional_centered)
{
	//DEPRECATED INTRO // REWRITE
	
	/*The purpose of this test is to check the run() function for the Solver1D class
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a Solver1D object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction Ez_)

	Lastly, the run() function is called.*/

	maxwell::Solver::Options solverOpts;
	
	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(TestMaxwellSolver1D::testModel, TestMaxwellSolver1D::defaultProbes, 
						TestMaxwellSolver1D::testSource, solverOpts);
	
	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
//
//TEST_F(TestMaxwellSolver1D, oneDimensional_upwind_PEC)
//{
//	int nx = 51;
//	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);
//
//	maxwell::Solver1D::Options solverOpts;
//
//	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
//
//	maxwell::Solver1D solver(solverOpts, mesh);
//	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
//	solver.setInitialField(FieldType::E, gaussianFunction);
//
//	Vector eOld = solver.getField(FieldType::E);
//	solver.run();
//	Vector eNew = solver.getField(FieldType::E);
//
//	double error = eOld.DistanceTo(eNew);
//	EXPECT_NEAR(0.0, error, 2e-3);
//
//}
//
//TEST_F(TestMaxwellSolver1D, oneDimensional_upwind_PMC)
//{
//	int nx = 51;
//	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);
//
//	maxwell::Solver1D::Options solverOpts;
//
//	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
//	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::PMC;
//
//	maxwell::Solver1D solver(solverOpts, mesh);
//	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
//	solver.setInitialField(FieldType::H, gaussianFunction);
//
//	Vector hOld = solver.getField(FieldType::H);
//	solver.run();
//	Vector hNew = solver.getField(FieldType::H);
//
//	double error = hOld.DistanceTo(hNew);
//	EXPECT_NEAR(0.0, error, 2e-3);
//
//}
//
//TEST_F(TestMaxwellSolver1D, oneDimensional_upwind_SMA)
//{
//	int nx = 51;
//	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);
//
//	maxwell::Solver1D::Options solverOpts;
//	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
//	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::SMA;
//	solverOpts.extractDataAtPoint = true;
//	IntegrationPoint ip;
//	double xPos = 1.0;
//	ip.Set1w(xPos, 0.0);
//	solverOpts.integPoint = ip;
//
//	maxwell::Solver1D solver(solverOpts, mesh);
//	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
//	solver.setInitialField(FieldType::E, gaussianFunction);
//
//	solver.run();
//
//	Vector eNew = solver.getField(FieldType::E);
//	Vector zero = eNew;
//	zero = 0.0;
//	double error = zero.DistanceTo(eNew);
//	EXPECT_NEAR(0.0, error, 2e-3);
//
//	Vector ePoint = solver.getFieldAtPoint();
//	for (int i = ePoint.Size() / 2; i < ePoint.Size(); i++) {
//		EXPECT_LE(ePoint[i],1.0);
//	}
//}
//
//TEST_F(TestMaxwellSolver1D, oneDimensional_two_materials)
//{
//	int nx = 100;
//	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);
//	HelperFunctions1D::setAttributeIntervalMesh(2, Vector({ 51,100 }), mesh);
//
//	maxwell::Solver1D::Options solverOpts;
//	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
//	solverOpts.t_final = 2.999;
//	
//	Vector eps = Vector({ 1.0, 2.0 });
//	Vector mu = Vector({ 1.0, 1.0 });
//	std::list<Material> matArray;
//	for (int i = 0; i < eps.Size(); i++) { //Check if mesh.attributes.Max() broken
//		Material matAux(eps[i], mu[i]);
//		matArray.push_back(matAux);
//	}
//
//	maxwell::Solver1D solver(solverOpts, mesh);
//	
//	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
//	solver.setInitialField(FieldType::E, gaussianFunctionHalfWidth);
//
//	Vector eOld = solver.getField(FieldType::E);
//	solver.run();
//	Vector eNew = solver.getField(FieldType::E);
//
//	double error = eOld.DistanceTo(eNew);
//	EXPECT_NEAR(0.0, error, 2e-3);
//}