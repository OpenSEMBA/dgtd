#include <gtest/gtest.h>

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"


using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class Solver2DTest : public ::testing::Test {
protected:
	static const int defaultNumberOfElements_X{ 3 };
	static const int defaultNumberOfElements_Y{ 3 };

	Model buildModel(
		const int nx = defaultNumberOfElements_X,
		const int ny = defaultNumberOfElements_Y,
		const Element::Type elType = Element::Type::TRIANGLE,
		const BdrCond& bdrB = BdrCond::PEC,
		const BdrCond& bdrR = BdrCond::PEC,
		const BdrCond& bdrT = BdrCond::PEC,
		const BdrCond& bdrL = BdrCond::PEC) {
		auto msh{ Mesh::MakeCartesian2D(nx,ny,elType) };
		return Model(msh, AttributeToMaterial{}, buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL));
	}

	Model buildModel(
		const int nx = defaultNumberOfElements_X,
		const int ny = defaultNumberOfElements_Y,
		const Element::Type elType = Element::Type::TRIANGLE,
		const double sx = 1.0,
		const double sy = 1.0,
		const BdrCond& bdrB = BdrCond::PEC,
		const BdrCond& bdrR = BdrCond::PEC,
		const BdrCond& bdrT = BdrCond::PEC,
		const BdrCond& bdrL = BdrCond::PEC) {
		auto msh{ Mesh::MakeCartesian2D(nx, ny, elType, false, sx, sy) };
		return Model(msh, AttributeToMaterial{}, buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL));
	}

	AttributeToBoundary buildAttrToBdrMap2D(const BdrCond& bdrB, const BdrCond& bdrR, const BdrCond& bdrT, const BdrCond& bdrL)
	{
		return {
			{1, bdrB},
			{2, bdrR},
			{3, bdrT},
			{4, bdrL},
		};
	}

	Vector fieldCenter{ { 0.5, 0.5 } };
};

class Solver2DSpectralTest : public Solver2DTest {

};

TEST_F(Solver2DTest, pec_centered_tris_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe(10)};
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	}; 

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//At the left boundary the electric field should be closed to zero and
	//the magnetic field reaches a maximum close to 1.0 or -1.0
	//(the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}

TEST_F(Solver2DSpectralTest, pec_centered_tris_1dot5D_spectral)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO(true)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}

TEST_F(Solver2DTest, pec_centered_quads_1dot5D)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DSpectralTest, pec_centered_quads_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}

TEST_F(Solver2DTest, pec_centered_quads_1dot5D_AMR)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.
	The field along Ez is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
	higher than the initial value.*/

	Mesh mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "amr-quad.mesh").c_str(),1,0) };
	mesh.UniformRefinement();

	AttributeToBoundary attToBdr{ {1, BdrCond::PMC}, {2,BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, Vector({0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 2e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DSpectralTest, pec_tris_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}

TEST_F(Solver2DTest, pec_quads_1dot5D_AMR)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.
	The field along Ez is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
	higher than the initial value.*/

	Mesh mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "amr-quad.mesh").c_str(),1,0)};
	mesh.UniformRefinement();

	AttributeToBoundary attToBdr{ {1, BdrCond::PMC}, {2,BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorConditions{} };

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, Vector({0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(2.5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 2e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DTest, pec_tris_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe(50) };

	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.4}},
		PointProbe{E, Z, {1.0, 0.4}},
		PointProbe{H, Y, {0.0, 0.4}},
		PointProbe{H, Y, {1.0, 0.4}}
	};

	maxwell::Solver solver{
	buildModel(
		9, 9, Element::Type::TRIANGLE, 2.0, 2.0,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC
	),
	probes,
	buildGaussianInitialField(E, 0.2, Vector({ 1.0,0.5 }), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-3)
		.setFinalTime(4.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DTest, pec_quads_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe(50) };

	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.4}},
		PointProbe{E, Z, {1.0, 0.4}},
		PointProbe{H, Y, {0.0, 0.4}},
		PointProbe{H, Y, {1.0, 0.4}}
	};

	maxwell::Solver solver{
	buildModel(
		10, 2, Element::Type::QUADRILATERAL, 1.0, 0.4, 
		BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC
	),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-3)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DSpectralTest, pec_quads_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}

TEST_F(Solver2DTest, sma_tris_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.4}},
		PointProbe{E, Z, {1.0, 0.4}},
		PointProbe{H, Y, {0.0, 0.4}},
		PointProbe{H, Y, {1.0, 0.4}}
	};


	maxwell::Solver solver{
		buildModel(
			10, 3, mfem::Element::Type::TRIANGLE, 1.0,3.0,
			BdrCond::PMC, BdrCond::SMA, BdrCond::PMC, BdrCond::SMA
		),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(2.0)
			.setOrder(3)
	};

	GridFunction eOld{ solver.getField(E,Z) };

	auto zeros{ eOld };
	zeros = 0.0;
	EXPECT_TRUE(eOld.DistanceTo(zeros) > 1e-2);

	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMin().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMin().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(2).findFrameWithMin().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DTest, sma_quads_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe(50) };
	//Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.4}},
		PointProbe{E, Z, {1.0, 0.4}},
		PointProbe{H, Y, {0.0, 0.4}},
		PointProbe{H, Y, {1.0, 0.4}}
	};

	maxwell::Solver solver{
	buildModel(
		11, 11, mfem::Element::Type::QUADRILATERAL,2.0, 2.0, 
		BdrCond::PMC, BdrCond::SMA, BdrCond::PMC, BdrCond::SMA
	),
	probes,
	buildGaussianInitialField(E, 0.15, Vector{{1.0,1.0}}, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-3)
		.setFinalTime(5.0)
		.setOrder(3)
	};

	GridFunction eOld{ solver.getField(E,Z) };

	auto zeros{ eOld };
	zeros = 0.0;
	EXPECT_TRUE( eOld.DistanceTo(zeros) > 1e-2);

	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMin().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMin().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(2).findFrameWithMin().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);

}

TEST_F(Solver2DTest, rotated_centered_quads_1dot5D)
{
	auto mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "severalrotatedquads.mesh").c_str(),1,0)};
	mesh.UniformRefinement();
	AttributeToBoundary attToBdr{ {1,BdrCond::PEC}, {2,BdrCond::PMC}};
	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorConditions{} };
	auto probes{ buildProbesWithAnExportProbe(30) };

	fieldCenter = Vector({ 2.0, 2.0 });

	maxwell::Solver solver {
		model,
		probes,
		buildGaussianInitialField(E, 0.5, fieldCenter, unitVec(Z), 1, Source::CartesianAngles({0.0,0.0,-M_PI_4})),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.95)
			.setCentered()
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver2DSpectralTest, periodic_centered_tris_spectral_and_base_comparison) {

	Probes probes;

	Mesh m;
	{
		Mesh square{ Mesh::MakeCartesian2D(9, 9, Element::TRIANGLE, false, 1.0, 1.0) };
		std::vector<Vector> translations{
			Vector({1.0, 0.0}),
			Vector({0.0, 1.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	maxwell::Solver solverSpectral{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	Vector zeroVec{ solver.getField(E,Z).Size() };
	zeroVec = 0.0;
	double tolerance{ 1e-5 };
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E,Z).Elem(i) - solverSpectral.getField(E,Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,Y).Elem(i) - solverSpectral.getField(H,Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,X).Elem(i) - solverSpectral.getField(H,X).Elem(i), tolerance);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E,Z).Elem(i) - solverSpectral.getField(E,Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,Y).Elem(i) - solverSpectral.getField(H,Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,X).Elem(i) - solverSpectral.getField(H,X).Elem(i), tolerance);
	}

}

TEST_F(Solver2DSpectralTest, periodic_tris_spectral_and_base_comparison) {

	Probes probes;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	maxwell::Solver solverSpectral{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	Vector zeroVec{ solver.getField(E,Z).Size() };
	zeroVec = 0.0;
	double tolerance{ 1e-5 };
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E,Z).Elem(i) - solverSpectral.getField(E,Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,Y).Elem(i) - solverSpectral.getField(H,Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,X).Elem(i) - solverSpectral.getField(H,X).Elem(i), tolerance);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E,Z).Elem(i) - solverSpectral.getField(E,Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,Y).Elem(i) - solverSpectral.getField(H,Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H,X).Elem(i) - solverSpectral.getField(H,X).Elem(i), tolerance);
	}

}

TEST_F(Solver2DTest, periodic_centered_tris)
{
	auto probes{ buildProbesWithAnExportProbe() };
	//Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};
	
	Mesh m;
	{
		Mesh square{ Mesh::MakeCartesian2D(9, 9, Element::TRIANGLE, false, 1.0, 1.0) };
		std::vector<Vector> translations{
			Vector({1.0, 0.0}),
			Vector({0.0, 1.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.1},
			Source::Position({0.5, 0.5}), // center_
			Source::Polarization(unitVec(Z)), // e polarization_
			Source::Propagation(unitVec(X))  // propagation direction
		),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-3 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave keeps travelling from left to right).
	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(2).findFrameWithMin().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(3).findFrameWithMin().second, tolerance);

}

TEST_F(Solver2DTest, periodic_centered_quads)
{

	auto probes{ buildProbesWithAnExportProbe(200) };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}
	};

	Mesh m;
	{
		Mesh square{ Mesh::MakeCartesian2D(9, 9, Element::QUADRILATERAL, false, 2.0, 2.0) };
		std::vector<Vector> translations{
			Vector({2.0, 0.0}),
			Vector({0.0, 2.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.2},
			Source::Position({1.0, 1.0}),           // center_
			Source::Polarization(unitVec(Z)),  // e polarization_
			Source::Propagation(unitVec(X))   // propagation direction
		),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 3e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave keeps travelling from left to right).
	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(2).findFrameWithMin().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(3).findFrameWithMin().second, tolerance);


}

TEST_F(Solver2DTest, periodic_tris)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	//Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {2.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {2.0, 0.5}}
	};

	Mesh m;
	{
		Mesh square{ Mesh::MakeCartesian2D(9, 9, Element::TRIANGLE, false, 2.0, 2.0) };
		std::vector<Vector> translations{
			Vector({2.0, 0.0}),
			Vector({0.0, 2.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.2},
			Source::Position({1.0, 1.0}), // center_
			Source::Polarization(unitVec(Z)), // e polarization_
			Source::Propagation(unitVec(X))  // propagation direction
		),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 3e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave keeps travelling from left to right).
	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(2).findFrameWithMin().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(3).findFrameWithMin().second, tolerance);

}

TEST_F(Solver2DTest, periodic_quads)
{
	auto probes{ buildProbesWithAnExportProbe(100) };
	//Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {2.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {2.0, 0.5}}
	};

	Mesh m;
	{
		Mesh square{ Mesh::MakeCartesian2D(9, 9, Element::QUADRILATERAL, false, 2.0, 2.0) };
		std::vector<Vector> translations{
			Vector({2.0, 0.0}),
			Vector({0.0, 2.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.2},
			Source::Position({1.0, 1.0}), // center_
			Source::Polarization(unitVec(Z)), // e polarization_
			Source::Propagation(unitVec(X))  // propagation direction
		),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 3e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave keeps travelling from left to right).
	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(2).findFrameWithMin().second, tolerance);
	EXPECT_NEAR(-1.0, solver.getPointProbe(3).findFrameWithMin().second, tolerance);

}

TEST_F(Solver2DTest, pec_centered_totalfieldinout_1dot5D)
{
	Mesh mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "4x4_Quadrilateral_InnerSquare_IntBdr.mesh").c_str(),1,0)};
	mesh.UniformRefinement();
	AttributeToBoundary attToBdr{ {1, BdrCond::PEC}, {2,BdrCond::PMC} };
	AttributeToInteriorConditions attToIntBdr{ {301,BdrCond::TotalFieldIn}, {302,BdrCond::TotalFieldOut} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe(20) };
	//probes.pointProbes = {
	//PointProbe{ E, Z, {0.5001, 0.5} },
	//PointProbe{ E, Z, {0.5, 0.5} },
	//PointProbe{ H, Y, {3.5, 0.5} },
	//PointProbe{ H, X, {3.5, 0.5} }
	//};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setCentered()
			.setFinalTime(4.0)
			.setOrder(2)
	};

	solver.run();

	//{
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().first, 1e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(3.5, solver.getPointProbe(2).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, 1e-3);
	//	EXPECT_NEAR(3.5, solver.getPointProbe(3).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, 1e-3);
	//}
}

TEST_F(Solver2DTest, sma_totalfieldinout_1dot5D)
{
	Mesh mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "4x4_Quadrilateral_1dot5D_IntBdr.mesh").c_str(),1,0)};
	AttributeToBoundary attToBdr{ {1, BdrCond::SMA}, {2,BdrCond::PMC} };
	AttributeToInteriorConditions attToIntBdr{ {301,BdrCond::TotalFieldIn}, {302,BdrCond::TotalFieldOut} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe(20) };
	//probes.pointProbes = {
	//PointProbe{ E, Z, {0.5001, 0.5} },
	//PointProbe{ E, Z, {0.5, 0.5} },
	//PointProbe{ H, Y, {3.5, 0.5} },
	//PointProbe{ H, X, {3.5, 0.5} }
	//};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(10.0)
			.setOrder(2)
	};

	solver.run();

	//{
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().first, 1e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(3.5, solver.getPointProbe(2).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, 1e-3);
	//	EXPECT_NEAR(3.5, solver.getPointProbe(3).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, 1e-3);
	//}
}

//TEST_F(Solver2DTest, DISABLED_quadraticMesh)
//{
//	Mesh mesh = Mesh::LoadFromFile("./testData/star-q2.mesh", 1, 0);
//	auto fec = std::make_unique<DG_FECollection>(4, 2, BasisType::GaussLobatto);
//	auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get());
//
//	Model model = Model(mesh, AttributeToMaterial{}, AttributeToBoundary{});
//
//	maxwell::Solver solver{
//		model,
//		buildProbesWithAnExportProbe(),
//		buildGaussianInitialField(E, Z, 0.4, mfem::Vector({-0.02566,0.03028})),
//		SolverOptions{}
//		.setTimeStep(5e-4)
//		.setFinalTime(2.0)
//		.setOrder(3)
//	};
//
//	auto normOld{ solver.getFields().getNorml2() };
//	solver.run();
//
//	double tolerance{ 1e-2 };
//	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);
//
//}
 
//TEST_F(Solver2DTest, DISABLED_periodic_strong) //TODO ADD ENERGY CHECK
//{
//	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
//	Vector periodic({ 0.0, 1.0 });
//	std::vector<Vector> trans;
//	trans.push_back(periodic);
//	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D, mesh2D.CreatePeriodicVertexMapping(trans));
//
//	maxwell::Solver::Options opts;
//	opts.evolution = FiniteElementEvolution::Options();
//	opts.evolution.disForm = DisForm::Strong;
//
//	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());
//
//	Probes probes;
//	probes.addExporterProbeToCollection(ExporterProbe());
//	probes.vis_steps = 20;
//
//	Sources sources;
//	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.0, 0.0 })));
//
//	maxwell::Solver solver(
//		model,
//		probes,
//		sources,
//		opts);
//
//	solver.run();
//
//}
//TEST_F(Solver2DTest, DISABLED_centered_NC_MESH) //TODO ADD ENERGY CHECK
//{
//	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
//	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.
//
//	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
//	A single 2D Gaussian on X and Y is declared along Ez.
//
//	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
//	with triangles and squares.
//	The field along Ez is extracted before and after the solver calls its run() method and evolves the
//	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
//	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
//	higher than the initial value.*/
//
//	const char* mesh_file = "star-mixed.mesh";
//	Mesh mesh(mesh_file);
//	mesh.UniformRefinement();
//	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());
//
//	Probes probes;
//	//probes.addExporterProbeToCollection(ExporterProbe());
//	//probes.vis_steps = 20;
//
//	Sources sources;
//	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 })));
//
//	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(2.92);
//	solverOpts.evolution.fluxType = FluxType::Centered;
//
//	maxwell::Solver solver(model, probes, sources, solverOpts);
//
//	GridFunction eOld = solver.getFieldInDirection(E, Z);
//	solver.run();
//	GridFunction eNew = solver.getFieldInDirection(E, Z);
//
//	EXPECT_GT(eOld.Max(), eNew.Max());
//}
//TEST_F(Solver2DTest, DISABLED_resonantBox)
//{
//	Mesh mesh2D = Mesh::MakeCartesian2D(21, 21, Element::Type::QUADRILATERAL);
//	std::vector<Attribute> attArrSingle = std::vector<Attribute>({ 1 });
//	Material mat11 = Material(1.0, 1.0);
//	std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
//	AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatMap(attArrSingle, matArrSimple);
//	AttributeToBoundary attToBdrVec;
//	Model model(mesh2D, attToMatVec, attToBdrVec);
//
//	double spread = 2.0;
//	double coeff = 20.0;
//	const Vector dev = Vector({ 0.0,0.0 });
//	Source EXFieldSource = Source(model, spread, coeff, dev, X, E); 
//	Sources sources;
//	sources.addSourceToVector(EXFieldSource);
//
//	Probes probes;
//	//probes.paraview = true;
//	probes.vis_steps = 100;
//
//	maxwell::Solver::Options solverOpts;
//
//	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
//	solverOpts.timeStep = 1e-4;
//	solverOpts.order = 1;
//
//	maxwell::Solver solver(model, probes,
//		sources, solverOpts);
//
//	solver.run();
//
//}

TEST_F(Solver2DTest, interiorPEC_sma_boundaries)
{
	Mesh mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "InteriorPEC2D.msh").c_str(),1,0)};
	mesh.UniformRefinement();
	AttributeToBoundary attToBdr{ {3,BdrCond::PMC}, {4,BdrCond::SMA } };
	AttributeToInteriorConditions attToIntBdr{ {2,BdrCond::PEC} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe(100) };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.2, Source::Position( {1.0, 0.0} ), unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver2DTest, normalStudy)
{
	Mesh mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "Maxwell2D_K2.mesh").c_str(),1,0)};
	AttributeToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model{ mesh, AttributeToMaterial{}, pecBdr, AttributeToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(100) };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.2, Source::Position({1.0, 0.0}), unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.0)
			.setOrder(1)
	};

	solver.run();

}

TEST_F(Solver2DTest, DISABLED_box_with_Gmsh)
{
	auto mesh = Mesh::LoadFromFile((gmshMeshesFolder() + "test.msh").c_str(), 1, 0);
	auto fec = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
	auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get(), 1, 0);
	auto model = Model(mesh, AttributeToMaterial{}, AttributeToBoundary{});

	auto probes{ buildProbesWithAnExportProbe(100) };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.1, mfem::Vector({0.5,0.5})),
		SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(4)
	};
}

TEST_F(Solver2DTest, 2D_box_resonant_mode)
{
	maxwell::Solver solver{ smbCase(getCaseName()) };

	solver.run();

	EXPECT_TRUE(false);
}