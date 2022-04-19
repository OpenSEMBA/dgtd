#include "gtest/gtest.h"
#include <math.h>

#include "maxwell1D/Solver.h"

using namespace maxwell1D;

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

}
using namespace AnalyticalFunctions1D;

class TestMaxwell1DSolver : public ::testing::Test {
protected:
	
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
};

TEST_F(TestMaxwell1DSolver, checkTwoAttributeMesh)
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
TEST_F(TestMaxwell1DSolver, oneDimensional_centered)
{
	/*The purpose of this test is to check the run() function for the Solver class
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a Solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction Ez_)

	Lastly, the run() function is called.*/

	int nx = 51;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);

	maxwell1D::Solver::Options solverOpts;
	
	solverOpts.evolutionOperatorOptions = FE_Evolution::Options();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell1D::Solver solver(solverOpts, mesh);
	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
	solver.setInitialField(FieldType::Electric, gaussianFunction);
	
	Vector eOld = solver.getField(FieldType::Electric);
	solver.run();
	Vector eNew = solver.getField(FieldType::Electric);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwell1DSolver, oneDimensional_upwind_PEC)
{
	int nx = 51;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);

	maxwell1D::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FE_Evolution::Options();

	maxwell1D::Solver solver(solverOpts, mesh);
	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
	solver.setInitialField(FieldType::Electric, gaussianFunction);

	Vector eOld = solver.getField(FieldType::Electric);
	solver.run();
	Vector eNew = solver.getField(FieldType::Electric);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwell1DSolver, oneDimensional_upwind_SMA)
{
	int nx = 51;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);

	maxwell1D::Solver::Options solverOpts;
	solverOpts.vis_steps = 25;
	solverOpts.paraview = true;

	solverOpts.evolutionOperatorOptions = FE_Evolution::Options();
	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::SMA;

	maxwell1D::Solver solver(solverOpts, mesh);
	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
	solver.setInitialField(FieldType::Electric, gaussianFunction);

	solver.run();
	Vector eNew = solver.getField(FieldType::Electric);
	Vector zero = eNew;
	zero = 0.0;

	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}