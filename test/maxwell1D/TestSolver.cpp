#include "gtest/gtest.h"
#include <math.h>

#include "maxwell1D/Solver.h"

using namespace Maxwell1D;

namespace AnalyticalFunctions1D {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;

	const double PI = atan(1.0) * 4;

	double gaussianFunction(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 2.0 * (pos[0] - center) /
			            (meshBoundingBoxMax[0] - meshBoundingBoxMin[0]);
		
		return exp(-20. * pow(normalizedPos, 2));
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

TEST_F(TestMaxwell1DSolver, oneDimensional)
{
	/*The purpose of this test is to check the run() function for the Solver class
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a Solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction Ez_)

	Lastly, the run() function is called.*/

	int nx = 11;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);

	Maxwell1D::Solver::Options opts;
	opts.order = 3;
	opts.dt = 1e-3;
	opts.t_final = 1000 * opts.dt;
	opts.vis_steps = 2;
	opts.paraview = true;

	Maxwell1D::Solver solver(opts, mesh);
	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);

	solver.setInitialElectricField(gaussianFunction);
	solver.run();
}

