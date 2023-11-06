#include <gtest/gtest.h>

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class Solver3DTest : public ::testing::Test {
protected:
	static const int defaultNumberOfElements_X{ 3 };
	static const int defaultNumberOfElements_Y{ 3 };
	static const int defaultNumberOfElements_Z{ 3 };

	Model buildModel(
		const int nx = defaultNumberOfElements_X,
		const int ny = defaultNumberOfElements_Y,
		const int nz = defaultNumberOfElements_Z,
		const Element::Type elType = Element::Type::HEXAHEDRON,
		const double sx = 1.0,
		const double sy = 1.0,
		const double sz = 1.0,
		
		const BdrCond& bdr1 = BdrCond::PEC,
		const BdrCond& bdr2 = BdrCond::PEC,
		const BdrCond& bdr3 = BdrCond::PEC,
		const BdrCond& bdr4 = BdrCond::PEC,
		const BdrCond& bdr5 = BdrCond::PEC,
		const BdrCond& bdr6 = BdrCond::PEC) {

		auto msh{ Mesh::MakeCartesian3D(nx, ny, nz, elType, sx, sy, sz) };

		return Model(msh, GeomTagToMaterial{}, buildAttrToBdrMap3D(bdr1, bdr2, bdr3, bdr4, bdr5, bdr6));
	}

	static GeomTagToBoundary buildAttrToBdrMap3D(const BdrCond& bdr1, const BdrCond& bdr2, const BdrCond& bdr3, const BdrCond& bdr4, const BdrCond& bdr5, const BdrCond& bdr6)
	{
		return {
			{1, bdr1},
			{2, bdr2},
			{3, bdr3},
			{4, bdr4},
			{5, bdr5},
			{6, bdr6},
		};
	}

	static void rotateMinus45degAlongXAxis(const Vector& oldP, Vector& newP)
	{
		assert(oldP.Size() == newP.Size());
		newP[0] = oldP[0];
		newP[1] = oldP[1] * cos(-M_PI / 4.0) - oldP[2] * sin(-M_PI / 4.0);
		newP[2] = oldP[1] * sin(-M_PI / 4.0) + oldP[2] * cos(-M_PI / 4.0);

	}

	static void rotateMinus45degAlongYAxis(const Vector& oldP, Vector& newP)
	{
		assert(oldP.Size() == newP.Size());
		newP[0] = oldP[0] *  cos(-M_PI / 4.0) + oldP[2] * sin(-M_PI / 4.0);
		newP[1] = oldP[1];
		newP[2] = oldP[0] * -sin(-M_PI / 4.0) + oldP[2] * cos(-M_PI / 4.0);

	}

	static void rotateMinus45degAlongZAxis(const Vector& oldP, Vector& newP)
	{
		assert(oldP.Size() == newP.Size());
		newP[0] = oldP[0] * cos(-M_PI / 4.0) - oldP[1] * sin(-M_PI / 4.0);
		newP[1] = oldP[0] * sin(-M_PI / 4.0) + oldP[1] * cos(-M_PI / 4.0);
		newP[2] = oldP[2];

	}
	static void rotateMinus90degAlongZAxis(const Vector& oldP, Vector& newP)
	{
		assert(oldP.Size() == newP.Size());
		newP[0] = oldP[0] * cos(-M_PI / 2.0) - oldP[1] * sin(-M_PI / 2.0);
		newP[1] = oldP[0] * sin(-M_PI / 2.0) + oldP[1] * cos(-M_PI / 2.0);
		newP[2] = oldP[2];
	}
};

class Solver3DSpectralTest : public Solver3DTest {

};

TEST_F(Solver3DTest, 3D_pec_centered_hexa_1dot5D)
{
	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
		buildModel(
			10,    1,   1, Element::Type::HEXAHEDRON, 
			1.0, 1.0, 1.0, 
			BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
			BdrCond::PMC,BdrCond::PEC,BdrCond::PEC
		),
		probes,
		buildGaussianInitialField(
			E, 0.1, 
			Source::Position({0.5,0.5,0.5}), 
			unitVec(Z)
		),
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

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DSpectralTest, 3D_pec_centered_hexa_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
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

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);


}

TEST_F(Solver3DSpectralTest, 3D_pec_centered_spectral_and_base_comparison)
{

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	Probes{},
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	maxwell::Solver solverSpectral{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	Probes{},
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	for (int i = 0; i < solver.getFields().allDOFs().Size(); ++i) {
		EXPECT_NEAR(solver.getFields().allDOFs()[i], solverSpectral.getFields().allDOFs()[i], 1e-5);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < solver.getFields().allDOFs().Size(); ++i) {
		EXPECT_NEAR(solver.getFields().allDOFs()[i], solverSpectral.getFields().allDOFs()[i], 1e-5);
	}

}

TEST_F(Solver3DTest, 3D_pec_upwind_hexa_1dot5D)
{
	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DSpectralTest, 3D_pec_upwind_hexa_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
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

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DTest, 3D_pec_centered_tetra_1dot5D)
{
	Probes probes{ buildProbesWithAnExportProbe(5) };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::TETRAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
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

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DSpectralTest, 3D_pec_centered_tetra_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
		buildModel(
			10, 1, 1, Element::Type::TETRAHEDRON,
			1.0, 1.0, 1.0,
			BdrCond::PEC, BdrCond::PMC, BdrCond::PEC,
			BdrCond::PMC, BdrCond::PEC, BdrCond::PEC),
			probes,
			buildGaussianInitialField(
				E, 0.1,
				Source::Position({ 0.5,0.5,0.5 }),
				unitVec(Z)
			),
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

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_pec_upwind_tetra_1dot5D)
{

	Probes probes{ buildProbesWithAnExportProbe(100) };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{E, Z, {1.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {1.0, 0.5, 0.5}}
	};

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::TETRAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5.0e-3)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DTest, 3D_gmsh_cube_upwind_tetra)
{
	Probes probes{ buildProbesWithAnExportProbe(20) };
	//probes.pointProbes = {
	//	PointProbe{E, Y, {0.0, 0.5, 0.5}},
	//	PointProbe{E, Y, {2.0, 0.5, 0.5}},
	//	PointProbe{H, Z, {0.0, 0.5, 0.5}},
	//	PointProbe{H, Z, {2.0, 0.5, 0.5}}
	//};
	
	Mesh m{ Mesh::LoadFromFile((gmshMeshesFolder() + "pureCube.msh").c_str(),1,0)};

	GeomTagToBoundary attToBdr{ {2,BdrCond::PMC},{3,BdrCond::PEC},{4,BdrCond::SMA} };
	Model model{ m, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{3.0}, 
			Source::Position    ({15.0, 0.5, 0.5}), 
			Source::Polarization(unitVec(Z)),
			Source::Propagation(unitVec(X)) 
		),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(15.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_pec_periodic_cube_centered_hexa)
{
	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Y, {0.0, 0.5, 0.5}},
		PointProbe{E, Y, {2.0, 0.5, 0.5}},
		PointProbe{H, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Z, {2.0, 0.5, 0.5}}
	};
	
	Mesh m;
	{
		Mesh cube{ Mesh::MakeCartesian3D(11,3,3,Element::HEXAHEDRON,2.0,2.0,2.0) };
		std::vector<Vector> translations{
			Vector({2.0, 0.0, 0.0}),
			Vector({0.0, 2.0, 0.0}),
			Vector({0.0, 0.0, 2.0})
		};
		m = Mesh::MakePeriodic(cube, cube.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.2}, 
			Source::Position    ({1.0, 0.5, 0.5}), 
			Source::Polarization(unitVec(Y)),
			Source::Propagation(unitVec(X)) 
		),
		SolverOptions{}
			.setTimeStep(7.5e-3)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_pec_periodic_cube_upwind_hexa)
{
	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Y, {0.0, 0.5, 0.5}},
		PointProbe{E, Y, {2.0, 0.5, 0.5}},
		PointProbe{H, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Z, {2.0, 0.5, 0.5}}
	};

	Mesh m;
	{
		Mesh cube{ Mesh::MakeCartesian3D(11,3,3,Element::HEXAHEDRON,2.0,2.0,2.0) };
		std::vector<Vector> translations{
			Vector({2.0, 0.0, 0.0}),
			Vector({0.0, 2.0, 0.0}),
			Vector({0.0, 0.0, 2.0})
		};
		m = Mesh::MakePeriodic(cube, cube.CreatePeriodicVertexMapping(translations));
	}

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.2},
			Source::Position({1.0, 0.5, 0.5}), // center_
			Source::Polarization(unitVec(Y)), // e polarization_
			Source::Propagation(unitVec(X)) 
		),
		SolverOptions{}
			.setTimeStep(3e-3)
			.setFinalTime(2.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DSpectralTest, 3D_pec_upwind_spectral_and_base_comparison) {

	maxwell::Solver solver{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	Probes{},
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	maxwell::Solver solverSpectral{
	buildModel(
		10,    1,   1, Element::Type::HEXAHEDRON,
		1.0, 1.0, 1.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	Probes{},
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	for (int i = 0; i < solver.getFields().allDOFs().Size(); ++i) {
		EXPECT_NEAR(solver.getFields().allDOFs()[i], solverSpectral.getFields().allDOFs()[i], 1e-5);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < solver.getFields().allDOFs().Size(); ++i) {
		EXPECT_NEAR(solver.getFields().allDOFs()[i], solverSpectral.getFields().allDOFs()[i], 1e-5);
	}
}

TEST_F(Solver3DTest, 3D_sma_upwind_hexa_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe(1) };
	//probes.pointProbes = {
	//	PointProbe{E, Z, {0.0, 0.5, 0.5}},
	//	PointProbe{E, Z, {1.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {0.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {1.0, 0.5, 0.5}}
	//};

	maxwell::Solver solver{
	buildModel(
		12,    4,   4, Element::Type::HEXAHEDRON,
		30.0, 15.0, 12.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::SMA,
		BdrCond::PMC,BdrCond::SMA,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 3.0,
		Source::Position({15.0,7.5,6.0}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5e-1)
		.setFinalTime(5.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DTest, feng_fss)
{
	auto probes{ buildProbesWithAnExportProbe(1000) };

	std::vector<double> pointR({ 0.01,-0.075,0.06 });
	std::vector<double> pointT({ 0.29,-0.075,0.06 });

	probes.pointProbes = {
		PointProbe{E, X, pointR},
		PointProbe{E, X, pointT},		
		PointProbe{E, Y, pointR},
		PointProbe{E, Y, pointT},		
		PointProbe{E, Z, pointR},
		PointProbe{E, Z, pointT},
		PointProbe{H, X, pointR},
		PointProbe{H, X, pointT},
		PointProbe{H, Y, pointR},
		PointProbe{H, Y, pointT},
		PointProbe{H, Z, pointR},
		PointProbe{H, Z, pointT}
	};

	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "fengfss.msh").c_str(),1,0)};
	mesh.Transform(rotateMinus90degAlongZAxis);
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	mfem::Vector center_(3);
	rotateMinus90degAlongZAxis(Vector({ 0.075,0.075,0.06 }), center_);
	mfem::Vector polarization_(3);
	rotateMinus90degAlongZAxis(unitVec(Z), polarization_);
	

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{0.015},
		Source::Position({ 0.075,0.075,0.06 }), // center_
		Source::Polarization(unitVec(Z)), // e polarization_
		Source::Propagation(unitVec(Y)) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(5e-7)
		.setFinalTime(0.50)
		.setOrder(1)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	for (int probeNumber = 0; probeNumber < probes.pointProbes.size(); probeNumber++) {
		std::ofstream file("tnf_" + std::to_string(probeNumber) + ".txt"); 
		file << "Time and " + std::to_string(probes.pointProbes[probeNumber].getFieldType()) + std::to_string(probes.pointProbes[probeNumber].getDirection()) + "\n";
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
			file << std::to_string(t) + " " + std::to_string(f) + "\n";
		}
	}

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, feng_fss_symmetry)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	std::vector<double> pointR({ 10.0, 37.5, 30 });
	std::vector<double> pointT({ 290.0, 37.5, 30 });

	probes.fieldProbes = {
		FieldProbe{pointR},
		FieldProbe{pointT}
	};

	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "Feng_FSS_Symmetry.msh").c_str(),1,0)};
	//mesh.Transform(rotateMinus90degAlongZAxis);
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA}};
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{16.0},
		Source::Position({ 75.0, 0.0, 0.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		Source::Propagation(unitVec(X)) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(1e-1)
		.setFinalTime(300.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	for (int probeNumber = 0; probeNumber < probes.fieldProbes.size(); probeNumber++) {
		std::ofstream file(getTestCaseName() + std::to_string(probeNumber) + ".txt");
		file << "Time // Ex // Ey // Ez // Hx // Hy // Hz //""\n";
		for (const auto& fm : solver.getFieldProbe(probeNumber).getFieldMovies()) {
			std::stringstream time, Ex, Ey, Ez, Hx, Hy, Hz;
			time << std::scientific << std::setprecision(7) << (fm.first); 
			Ex << std::scientific << std::setprecision(7) << fm.second.Ex; Ey << std::scientific << std::setprecision(7) << fm.second.Ey; Ez << std::scientific << std::setprecision(7) << fm.second.Ez; 
			Hx << std::scientific << std::setprecision(7) << fm.second.Hx; Hy << std::scientific << std::setprecision(7) << fm.second.Hy; Hz << std::scientific << std::setprecision(7) << fm.second.Hz;
			file << time.str() + " " + Ex.str() + " " + Ey.str() + " " + Ez.str() + " " + Hx.str() + " " + Hy.str() + " " + Hz.str() + "\n";
		}
	}

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, feng_fss_manual)
{
	auto mesh{ Mesh::LoadFromFile((mfemMeshes3DFolder() + "fengfssmanual.mesh").c_str(),1,0)};
	Array<Refinement> refinement_list;
	refinement_list.Append(Refinement(0, 2));
	refinement_list.Append(Refinement(1, 2));
	refinement_list.Append(Refinement(2, 2));
	refinement_list.Append(Refinement(3, 2));
	refinement_list.Append(Refinement(4, 2));
	refinement_list.Append(Refinement(5, 2));
	refinement_list.Append(Refinement(6, 2));
	refinement_list.Append(Refinement(7, 2));
	mesh.GeneralRefinement(refinement_list);
	refinement_list.Append(Refinement(8, 2));
	refinement_list.Append(Refinement(9, 2));
	refinement_list.Append(Refinement(10, 2));
	refinement_list.Append(Refinement(11, 2));
	refinement_list.Append(Refinement(12, 2));
	refinement_list.Append(Refinement(13, 2));
	refinement_list.Append(Refinement(14, 2));
	refinement_list.Append(Refinement(15, 2));
	mesh.GeneralRefinement(refinement_list);
	//refinement_list.Append(Refinement(16, 2));
	//refinement_list.Append(Refinement(17, 2));
	//refinement_list.Append(Refinement(18, 2));
	//refinement_list.Append(Refinement(19, 2));
	mesh.GeneralRefinement(refinement_list);
	mesh.Transform(rotateMinus90degAlongZAxis);
	GeomTagToBoundary attToBdr{ 
		{2, BdrCond::PEC},
		{3, BdrCond::PMC},
		{4, BdrCond::SMA}
	};
	GeomTagToInteriorConditions attToIntBdr{ {5, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, attToIntBdr };

	mfem::Vector center(3);
	rotateMinus90degAlongZAxis(Vector({ 0.15,0.15,0.06 }), center);
	mfem::Vector polarization(3);
	rotateMinus90degAlongZAxis(unitVec(Z), polarization);

	auto probes{ buildProbesWithAnExportProbe(1000) };

	std::vector<double> pointR({ 0.01,-0.075,0.06 });
	std::vector<double> pointT({ 0.29,-0.075,0.06 });

	probes.pointProbes = {
		PointProbe{E, X, pointR},
		PointProbe{E, X, pointT},
		PointProbe{E, Y, pointR},
		PointProbe{E, Y, pointT},
		PointProbe{E, Z, pointR},
		PointProbe{E, Z, pointT},
		PointProbe{H, X, pointR},
		PointProbe{H, X, pointT},
		PointProbe{H, Y, pointR},
		PointProbe{H, Y, pointT},
		PointProbe{H, Z, pointR},
		PointProbe{H, Z, pointT}
	};

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{0.015},
		Source::Position({ 0.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		Source::Propagation(unitVec(X)) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(1e-7)
		.setFinalTime(0.0001)
		.setOrder(1)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	
	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, interiorPEC_sma_boundaries)
{
	Mesh mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "InteriorPEC3D.msh").c_str(),1,0)};
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC}, {4,BdrCond::SMA } };
	GeomTagToInteriorConditions attToIntBdr{ {5,BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{7.5},
			Source::Position({ 0.0 }), // center
			Source::Polarization(unitVec(Z)), // e polarization
			mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
	),
		SolverOptions{}
			.setTimeStep(1.0)
			.setFinalTime(30.0)
			.setOrder(2)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, interiorPEC_fss_hexas)
{
	auto probes{ buildProbesWithAnExportProbe(2) };

	std::vector<double> pointR({ 25, 25, 25 });
	std::vector<double> pointT({ 275, 25, 25 });

	probes.pointProbes = {
		PointProbe{E, X, pointR},
		PointProbe{E, X, pointT},
		PointProbe{E, Y, pointR},
		PointProbe{E, Y, pointT},
		PointProbe{E, Z, pointR},
		PointProbe{E, Z, pointT},
		PointProbe{H, X, pointR},
		PointProbe{H, X, pointT},
		PointProbe{H, Y, pointR},
		PointProbe{H, Y, pointT},
		PointProbe{H, Z, pointR},
		PointProbe{H, Z, pointT}
	};

	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "fsshexas.msh").c_str(),1,0)};
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{16},
		Source::Position({ 70, 0.0, 0.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector(unitVec(X)) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(7.5e-1)
		.setFinalTime(270.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	for (int probeNumber = 0; probeNumber < probes.pointProbes.size(); probeNumber++) {
		std::ofstream file("fss_sym_" + std::to_string(probeNumber) + ".txt");
		file << "Time and " + std::to_string(probes.pointProbes[probeNumber].getFieldType()) + std::to_string(probes.pointProbes[probeNumber].getDirection()) + "\n";
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
			file << std::to_string(t) + " " + std::to_string(f) + "\n";
		}
	}

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, 3D_minimal_tetra)
{
	auto probes{ buildProbesWithAnExportProbe() };
	Mesh mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "twotetras.msh").c_str(),1,0)};

	GeomTagToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh, GeomTagToMaterial{}, GeomTagToBoundary{}, GeomTagToInteriorConditions{} };

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{0.16},
		Source::Position({ 0.7, 0.0, 0.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector(unitVec(X)) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, sma_upwind_hex_1dot5D)
{

	Probes probes{ buildProbesWithAnExportProbe(200) };
	//probes.pointProbes = {
	//	PointProbe{E, Z, {0.0, 0.5, 0.5}},
	//	PointProbe{E, Z, {1.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {0.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {1.0, 0.5, 0.5}}
	//};

	maxwell::Solver solver{
		buildModel(
			10, 2, 2, Element::Type::HEXAHEDRON,
			1.0, 0.4, 0.4,
			BdrCond::PEC, BdrCond::PMC, BdrCond::SMA,
			BdrCond::PMC, BdrCond::SMA, BdrCond::PEC),
			probes,
			buildGaussianInitialField(
				E, 0.1,
				Source::Position({ 0.5,0.5,0.5 }),
				unitVec(Z)
			),
			SolverOptions{}
			.setTimeStep(1e-4)
			.setFinalTime(2.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DTest, 3D_pec_centered_hexa_totalfieldin)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {3.0, 0.5, 0.5}},
		PointProbe{E, Z, {5.0, 0.5, 0.5}},
		PointProbe{H, Y, {3.0, 0.5, 0.5}},
		PointProbe{H, Y, {5.0, 0.5, 0.5}}
	};
	auto mesh{ Mesh::LoadFromFile((mfemMeshes3DFolder() + "beam_hex_totalfieldin.mesh").c_str(), 1, 0) };
	GeomTagToBoundary att2bdr{ {1, BdrCond::PMC}, {2, BdrCond::PEC}, {3, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, GeomTagToInteriorConditions());

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(8.0)
			.setOrder(3)
	};

	solver.run();
	
	{
		auto frame{ solver.getPointProbe(0).getFieldMovie()};
		auto expected_t = 6.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(1).getFieldMovie() };
		auto expected_t = 8.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(2).getFieldMovie() };
		auto expected_t = 6.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(3).getFieldMovie() };
		auto expected_t = 8.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}

}

TEST_F(Solver3DTest, feng_fss_flat)
{
	auto probes{ buildProbesWithAnExportProbe(50) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_Feng_FSS_Flat.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary attToBdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::SMA} };
	GeomTagToInteriorConditions att2IntCond{ {60, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, att2IntCond };

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianPlanewave(0.010, 0.1, unitVec(Y), unitVec(X)),
	SolverOptions{}
		.setTimeStep(9e-4)
		.setFinalTime(1.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, 3D_pec_upwind_box_totalfieldscatteredfield)
{
	auto probes{ buildProbesWithAnExportProbe(30) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_MinimalistBox.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::SMA} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, GeomTagToInteriorConditions());

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.5, 5.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(10.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, 3D_pec_centered_beam_totalfieldscatteredfield)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {2.0, 0.5, 0.5}},
		PointProbe{E, Z, {5.0, 0.5, 0.5}},
		PointProbe{H, Y, {2.0, 0.5, 0.5}},
		PointProbe{H, Y, {5.0, 0.5, 0.5}}
	};
	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PMC}, {1, BdrCond::PEC}, {3, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, GeomTagToInteriorConditions());

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, upwind_beam_totalfieldscatteredfield_inout)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {2.0, 0.5, 0.5}},
		PointProbe{E, Z, {5.0, 0.5, 0.5}},
		PointProbe{H, Y, {2.0, 0.5, 0.5}},
		PointProbe{H, Y, {5.0, 0.5, 0.5}}
	};
	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PMC}, {1, BdrCond::PEC}, {3, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, GeomTagToInteriorConditions());

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, dualintbdr_upwind_beam_totalfieldscatteredfield_in)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_DualSurface_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC}, {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, GeomTagToInteriorConditions());

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 2.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(10.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, 3D_pec_centered_innerbox_totalfieldinout)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_Box.msh").c_str(), 1, 0, true) };
	mesh.UniformRefinement();
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, GeomTagToInteriorConditions());

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(10.0)
			.setOrder(3)
	};

	solver.run();

}


TEST_F(Solver3DTest, centered_beam_totalfieldscatteredfield_inout_intbdr)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, att2IntCond);

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, centered_beam_totalfieldscatteredfield_inout_intbdr_RtL)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam_RtL.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, att2IntCond);

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 8.0, unitVec(Z), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, upwind_beam_totalfieldscatteredfield_inout_intbdr)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, att2IntCond);

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, upwind_beam_totalfieldscatteredfield_inout_intbdr_RtL)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam_RtL.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), att2bdr, att2IntCond);

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 8.0, unitVec(Z), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}