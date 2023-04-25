#include "gtest/gtest.h"

#include "AnalyticalFunctions3D.h"
#include "SourceFixtures.h"
#include "maxwell/Solver.h"

#include <iostream>
#include <fstream>

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;
using namespace AnalyticalFunctions3D;

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

		return Model(msh, AttributeToMaterial{}, buildAttrToBdrMap3D(bdr1, bdr2, bdr3, bdr4, bdr5, bdr6));
	}

	static AttributeToBoundary buildAttrToBdrMap3D(const BdrCond& bdr1, const BdrCond& bdr2, const BdrCond& bdr3, const BdrCond& bdr4, const BdrCond& bdr5, const BdrCond& bdr6)
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

	static Probes buildProbesWithAnExportProbe()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
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

TEST_F(Solver3DTest, 3D_pec_centered_hexa_1dot5D_spectral)
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

TEST_F(Solver3DTest, 3D_pec_centered_spectral_and_base_comparison) 
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

	for (int i = 0; i < solver.getFields().E[X].Size(); ++i) {
		EXPECT_NEAR(solver.getFields().E[X].Elem(i), solverSpectral.getFields().E[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Y].Elem(i), solverSpectral.getFields().E[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Z].Elem(i), solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[Y].Elem(i), solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < solver.getFields().E[X].Size(); ++i) {
		EXPECT_NEAR(solver.getFields().E[X].Elem(i), solverSpectral.getFields().E[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Y].Elem(i), solverSpectral.getFields().E[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Z].Elem(i), solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[Y].Elem(i), solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
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

TEST_F(Solver3DTest, 3D_pec_upwind_hexa_1dot5D_spectral)
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
	Probes probes;
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

TEST_F(Solver3DTest, 3D_pec_centered_tetra_1dot5D_spectral)
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

TEST_F(Solver3DTest, sma_upwind_hex_1dot5D)
{

	Probes probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 200;
	//probes.pointProbes = {
	//	PointProbe{E, Z, {0.0, 0.5, 0.5}},
	//	PointProbe{E, Z, {1.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {0.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {1.0, 0.5, 0.5}}
	//};

	maxwell::Solver solver{
	buildModel(
		10,    2,   2, Element::Type::HEXAHEDRON,
		1.0, 0.4, 0.4,
		BdrCond::PEC,BdrCond::PMC,BdrCond::SMA,
		BdrCond::PMC,BdrCond::SMA,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.1,
		Source::Position({0.5,0.5,0.5}),
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

TEST_F(Solver3DTest, 3D_pec_upwind_tetra_1dot5D)
{

	Probes probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 5;
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
			E,
			Source::Position    ({1.0, 0.5, 0.5}), // center
			Source::Polarization({0.0, 1.0, 0.0}), // e polarization
			mfem::Vector        ({1.0, 0.0, 0.0})  // propagation direction
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
			E,
			Source::Position({1.0, 0.5, 0.5}), // center
			Source::Polarization({0.0, 1.0, 0.0}), // e polarization
			mfem::Vector({1.0, 0.0, 0.0})  // propagation direction
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

TEST_F(Solver3DTest, 3D_rotated_M45X_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(
	10,    1,   1, Element::HEXAHEDRON,
	1.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongXAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	Probes probes;
	probes.pointProbes = {
	PointProbe{H, Y, {0.0, 0.5, 0.5}},
	PointProbe{H, Y, {1.0, 0.5, 0.5}},
	PointProbe{H, Z, {0.0, 0.5, 0.5}},
	PointProbe{H, Z, {1.0, 0.5, 0.5}},
	PointProbe{E, Y, {0.0, 0.5, 0.5}},
	PointProbe{E, Y, {1.0, 0.5, 0.5}},
	PointProbe{E, Z, {0.0, 0.5, 0.5}},
	PointProbe{E, Z, {1.0, 0.5, 0.5}},
	};

	mfem::Vector center(3);
	rotateMinus45degAlongXAxis(Vector({ 0.5,0.5,0.5 }), center);
	mfem::Vector polarization(3);
	rotateMinus45degAlongXAxis(unitVec(Z), polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, center, polarization, 1, Source::CartesianAngles({ M_PI_4, 0.0, 0.0 })),
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
	EXPECT_NEAR(0.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(sqrt(2.0) / 2.0, solver.getPointProbe(4).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(sqrt(2.0) / 2.0, solver.getPointProbe(5).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(sqrt(2.0) / 2.0, solver.getPointProbe(6).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(sqrt(2.0) / 2.0, solver.getPointProbe(7).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_rotated_M45Y_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(
	10,    1,   1, Element::HEXAHEDRON,
	1.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongYAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorBoundary{} };

	Probes probes;
	probes.pointProbes = {
	PointProbe{E, X, {-0.35355339059327, 0.5, 0.35355339059327}},
	PointProbe{E, X, { 0.35355339059327, 0.5, 1.0606601717798 }},
	PointProbe{E, Z, {-0.35355339059327, 0.5, 0.35355339059327}},
	PointProbe{E, Z, { 0.35355339059327, 0.5, 1.0606601717798 }}
	};

	mfem::Vector center(3);
	rotateMinus45degAlongYAxis(Vector({ 0.5,0.5,0.5 }), center);
	mfem::Vector polarization(3);
	rotateMinus45degAlongYAxis(unitVec(Z), polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, center, polarization, 1, Source::CartesianAngles({ 0.0, M_PI_4, 0.0 })),
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

	EXPECT_NEAR(-sqrt(2.0) / 2.0, solver.getPointProbe(0).findFrameWithMin().second, tolerance);
	EXPECT_NEAR(-sqrt(2.0) / 2.0, solver.getPointProbe(1).findFrameWithMin().second, tolerance);
	EXPECT_NEAR(sqrt(2.0) / 2.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(sqrt(2.0) / 2.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_rotated_M45Z_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(
	10,    1,   1, Element::HEXAHEDRON,
	1.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongZAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorBoundary{} };

	Probes probes;
	probes.pointProbes = {
	PointProbe{E, Z, {0.35355339059327, 0.35355339059327, 0.5}},
	PointProbe{E, Z, {1.0606601717798, -0.35355339059327, 0.5}},
	PointProbe{H, Y, {0.35355339059327, 0.35355339059327, 0.5}},
	PointProbe{H, Y, {1.0606601717798, -0.35355339059327, 0.5}}
	};

	mfem::Vector center(3);
	rotateMinus45degAlongZAxis(Vector({ 0.5,0.5,0.5 }), center);
	mfem::Vector polarization(3);
	rotateMinus45degAlongZAxis(unitVec(Z), polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, center, polarization, 1, Source::CartesianAngles({ 0.0, 0.0, M_PI_4 })),
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

	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(1.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_rotated_AllDir_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(
		10,    1,   1, Element::HEXAHEDRON, 
		1.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongXAxis);
	mesh.Transform(rotateMinus45degAlongYAxis);
	mesh.Transform(rotateMinus45degAlongZAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	Probes probes;
	probes.pointProbes = {
	PointProbe{E, X, {0.5, 0.5, 0.0}},
	PointProbe{E, X, {1.0, 0.0, 0.707}},
	PointProbe{E, Y, {0.5, 0.5, 0.0}},
	PointProbe{E, Y, {1.0, 0.0, 0.707}},
	PointProbe{E, Z, {0.5, 0.5, 0.0}},
	PointProbe{E, Z, {1.0, 0.0, 0.707}},
	};


	mfem::Vector center({ 0.5,0.5,0.5 }), rCenter1(3), rCenter2(3), rotCenter(3);
	mfem::Vector polarization(3), tPolarization1(3), tPolarization2(3);

	rotateMinus45degAlongXAxis(center, rCenter1);
	rotateMinus45degAlongYAxis(rCenter1, rCenter2);
	rotateMinus45degAlongZAxis(rCenter2, rotCenter);

	rotateMinus45degAlongXAxis(unitVec(Z), tPolarization1);
	rotateMinus45degAlongYAxis(tPolarization1, tPolarization2);
	rotateMinus45degAlongZAxis(tPolarization2, polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, rotCenter, polarization, 1, Source::CartesianAngles({ M_PI_4, M_PI_4, M_PI_4 })),
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

	EXPECT_NEAR(0.146486, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.146486, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.853781, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.853781, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.5     , solver.getPointProbe(4).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.5     , solver.getPointProbe(5).findFrameWithMax().second, tolerance);

}

TEST_F(Solver3DTest, 3D_rotated_AllDir_upwind_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(
		10,    1,   1, Element::HEXAHEDRON,
		1.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongXAxis);
	mesh.Transform(rotateMinus45degAlongYAxis);
	mesh.Transform(rotateMinus45degAlongZAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	Probes probes;
	probes.pointProbes = {
	PointProbe{E, X, {0.5, 0.5, 0.0}},
	PointProbe{E, X, {1.0, 0.0, 0.707}},
	PointProbe{E, Y, {0.5, 0.5, 0.0}},
	PointProbe{E, Y, {1.0, 0.0, 0.707}},
	PointProbe{E, Z, {0.5, 0.5, 0.0}},
	PointProbe{E, Z, {1.0, 0.0, 0.707}},
	};


	mfem::Vector center({ 0.5,0.5,0.5 }), rCenter1(3), rCenter2(3), rotCenter(3);
	mfem::Vector polarization(3), tPolarization1(3), tPolarization2(3);

	rotateMinus45degAlongXAxis(center, rCenter1);
	rotateMinus45degAlongYAxis(rCenter1, rCenter2);
	rotateMinus45degAlongZAxis(rCenter2, rotCenter);

	rotateMinus45degAlongXAxis(unitVec(Z), tPolarization1);
	rotateMinus45degAlongYAxis(tPolarization1, tPolarization2);
	rotateMinus45degAlongZAxis(tPolarization2, polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, rotCenter, polarization, 1, Source::CartesianAngles({ M_PI_4, M_PI_4, M_PI_4 })),
	SolverOptions{}
		.setTimeStep(5e-3)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	EXPECT_NEAR(0.146486, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.146486, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.853781, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.853781, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.5, solver.getPointProbe(4).findFrameWithMax().second, tolerance);
	EXPECT_NEAR(0.5, solver.getPointProbe(5).findFrameWithMax().second, tolerance);

}
TEST_F(Solver3DTest, 3D_pec_upwind_spectral_and_base_comparison) {

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

	for (int i = 0; i < solver.getFields().E[X].Size(); ++i) {
		EXPECT_NEAR(solver.getFields().E[X].Elem(i), solverSpectral.getFields().E[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Y].Elem(i), solverSpectral.getFields().E[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Z].Elem(i), solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[Y].Elem(i), solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < solver.getFields().E[X].Size(); ++i) {
		EXPECT_NEAR(solver.getFields().E[X].Elem(i), solverSpectral.getFields().E[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Y].Elem(i), solverSpectral.getFields().E[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().E[Z].Elem(i), solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[Y].Elem(i), solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(solver.getFields().H[X].Elem(i), solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}
}

TEST_F(Solver3DTest, 3D_sma_upwind_hexa_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 1;
	//probes.pointProbes = {
	//	PointProbe{E, Z, {0.0, 0.5, 0.5}},
	//	PointProbe{E, Z, {1.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {0.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {1.0, 0.5, 0.5}}
	//};

	maxwell::Solver solver{
	buildModel(
		8,    4,   4, Element::Type::HEXAHEDRON,
		0.3, 0.15, 0.12,
		BdrCond::PEC,BdrCond::PMC,BdrCond::SMA,
		BdrCond::PMC,BdrCond::SMA,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.03,
		Source::Position({0.15,0.075,0.06}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(1e-5)
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


TEST_F(Solver3DTest, feng_fss)
{
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 1000;

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

	auto mesh{ Mesh::LoadFromFile("./testData/fengfss.msh",1,0) };
	mesh.Transform(rotateMinus90degAlongZAxis);
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorBoundary{} };

	mfem::Vector center(3);
	rotateMinus90degAlongZAxis(Vector({ 0.075,0.075,0.06 }), center);
	mfem::Vector polarization(3);
	rotateMinus90degAlongZAxis(unitVec(Z), polarization);
	

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{0.015},
		E,
		Source::Position({ 0.075,0.075,0.06 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
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
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 1000;

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

	auto mesh{ Mesh::LoadFromFile("./TestData/fengfssflatsym.msh",1,0) };
	mesh.Transform(rotateMinus90degAlongZAxis);
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA}};
	AttributeToInteriorBoundary attToIntBdr{ {5, BdrCond::PEC }, {6, BdrCond::NONE} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr};

	mfem::Vector center(3);
	rotateMinus90degAlongZAxis(Vector({ 0.075,0.075,0.06 }), center);
	mfem::Vector polarization(3);
	rotateMinus90degAlongZAxis(unitVec(Z), polarization);


	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{0.015},
		E,
		Source::Position({ 0.075,0.075,0.06 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(8e-7)
		.setFinalTime(0.05)
		.setOrder(1)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	for (int probeNumber = 0; probeNumber < probes.pointProbes.size(); probeNumber++) {
		std::ofstream file("tnf_sym_" + std::to_string(probeNumber) + ".txt");
		file << "Time and " + std::to_string(probes.pointProbes[probeNumber].getFieldType()) + std::to_string(probes.pointProbes[probeNumber].getDirection()) + "\n";
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
			file << std::to_string(t) + " " + std::to_string(f) + "\n";
		}
	}

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, feng_fss_manual)
{
	auto mesh{ Mesh::LoadFromFile("./testData/fengfssmanual.mesh",1,0) };
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
	AttributeToBoundary attToBdr{ 
		{2, BdrCond::PEC},
		{3, BdrCond::PMC},
		{4, BdrCond::SMA}
	};
	AttributeToInteriorBoundary attToIntBdr{ {5, BdrCond::PEC} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	mfem::Vector center(3);
	rotateMinus90degAlongZAxis(Vector({ 0.15,0.15,0.06 }), center);
	mfem::Vector polarization(3);
	rotateMinus90degAlongZAxis(unitVec(Z), polarization);

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 1000;

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
		E,
		Source::Position({ 0.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
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
	Mesh mesh{ Mesh::LoadFromFile("./testData/InteriorPEC3D.msh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC}, {4,BdrCond::SMA } };
	AttributeToInteriorBoundary attToIntBdr{ {5,BdrCond::PEC} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 100;

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.15},
			E,
			Source::Position({ 0.0 }), // center
			Source::Polarization(unitVec(Z)), // e polarization
			mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
	),
		SolverOptions{}
			.setTimeStep(2.5e-4)
			.setFinalTime(1.0)
			.setOrder(2)
	};

	solver.run();

}

TEST_F(Solver3DTest, interiorPEC_fss_hexas)
{
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 500;

	auto mesh{ Mesh::LoadFromFile("./TestData/fsshexasmoredetail.msh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA} };
	AttributeToInteriorBoundary attToIntBdr{ {5, BdrCond::PEC } };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{1.6},
		E,
		Source::Position({ 9.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(1e-4)
		.setFinalTime(15.0)
		.setOrder(2)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}
