#include "gtest/gtest.h"

#include "maxwell/Solver.h"

using namespace maxwell;
using namespace mfem;

class TestSolver3D : public ::testing::Test {
protected:

};

TEST_F(TestSolver3D, DISABLED_threeDimensional)
{
	Mesh mesh = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::HEXAHEDRON);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Sources sources{
		GaussianInitialField{E, Z, 0.2, 200.0, Vector({ 0.0,0.0,0.0 })}
	};

	Probes probes;

	maxwell::Solver solver(
		model,
		probes,
		sources,
		SolverOptions().setFinalTime(0.5));

	solver.run();
}
