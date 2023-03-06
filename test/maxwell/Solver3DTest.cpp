#include "gtest/gtest.h"

#include "AnalyticalFunctions3D.h"
#include "SourceFixtures.h"
#include "maxwell/Solver.h"

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

		return Model(Mesh::MakeCartesian3D(nx, ny, nz, elType, sx, sy, sz), AttributeToMaterial{}, buildAttrToBdrMap3D(bdr1, bdr2, bdr3, bdr4, bdr5, bdr6));
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
};

TEST_F(Solver3DTest, centered_hexa_1dot5D)
{
	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(
		10,1,1, Element::Type::HEXAHEDRON, 
		3.0, 1.0, 1.0, 
		BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
		BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(
		E, 0.2, 
		Source::Position({1.5,0.5,0.5}), 
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, centered_hexa_1dot5D_spectral)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
		.setSpectralEO()
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, upwind_hexa_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver3DTest, upwind_hexa_1dot5D_spectral)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(6.0)
		.setOrder(3)
		.setSpectralEO()
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver3DTest, centered_tetra_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 100;

	maxwell::Solver solver{
	buildModel(15,1,1, Element::Type::TETRAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-3)
		.setCentered()
		.setFinalTime(1.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, centered_tetra_1dot5D_spectral)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 500;

	maxwell::Solver solver{
	buildModel(15,1,1, Element::Type::TETRAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-4)
		.setCentered()
		.setFinalTime(1.0)
		.setOrder(3)
		.setSpectralEO(false, 0, true)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, upwind_tetra_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 500;

	maxwell::Solver solver{
	buildModel(10,1,1, Element::Type::TETRAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(30.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver3DTest, periodic_x_centered_tetra_1dot5)
{
	Mesh m{ Mesh::MakeCartesian3D(11,3,3,Element::TETRAHEDRON, 3.0, 3.0, 3.0) };

	Vector xTr({ 3.0,3.0,3.0 });
	std::vector<Vector> translations{ xTr };

	Mesh mPer{ Mesh::MakePeriodic(m, m.CreatePeriodicVertexMapping(translations)) };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	AttributeToBoundary attToBdr{ { 1,BdrCond::PEC }, { 2,BdrCond::PMC }, { 3,BdrCond::PMC }, { 4,BdrCond::PMC }, { 5,BdrCond::PMC }, { 6,BdrCond::PEC } };

	Model model{ mPer,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.1, mfem::Vector({0.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(1e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(2)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, periodic_cube_centered)
{
	Mesh m;
	{
		Mesh cube{ Mesh::MakeCartesian3D(10,3,3,Element::HEXAHEDRON,4.0,1.0,1.0) };
		std::vector<Vector> translations{
			Vector({4.0, 0.0, 0.0}),
			Vector({0.0, 1.0, 0.0}),
			Vector({0.0, 0.0, 1.0})
		};
		m = Mesh::MakePeriodic(cube, cube.CreatePeriodicVertexMapping(translations));
	}
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 100;

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.25}, 
			Source::Position    ({1.0, 0.5, 0.5}), // center
			Source::Polarization({0.0, 1.0, 0.0}), // e polarization
			mfem::Vector        ({1.0, 0.0, 0.0})  // propagation direction
		),
		SolverOptions{}
			.setTimeStep(5e-4)
			.setCentered()
			.setFinalTime(4.0)
			.setOrder(3)
	};

	solver.run();
}

TEST_F(Solver3DTest, periodic_cube_upwind)
{
	Mesh m;
	{
		Mesh cube{ Mesh::MakeCartesian3D(10,3,3,Element::HEXAHEDRON,4.0,1.0,1.0) };
		std::vector<Vector> translations{
			Vector({4.0, 0.0, 0.0}),
			Vector({0.0, 1.0, 0.0}),
			Vector({0.0, 0.0, 1.0})
		};
		m = Mesh::MakePeriodic(cube, cube.CreatePeriodicVertexMapping(translations));
	}
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 100;

	Model model{ m };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.25},
			Source::Position({1.0, 0.5, 0.5}), // center
			Source::Polarization({0.0, 1.0, 0.0}), // e polarization
			mfem::Vector({1.0, 0.0, 0.0})  // propagation direction
		),
		SolverOptions{}
			.setTimeStep(5e-4)
			.setFinalTime(4.0)
			.setOrder(3)
	};

	solver.run();
}

TEST_F(Solver3DTest, DISABLED_box_pec_upwind)
{
	/*The purpose of this test is to check the run() function for the solver object
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction ez_)

	Lastly, the run() function is called.*/

	//maxwell::Solver solver{
	//buildStandardModel(3,3,3),
	//buildProbesWithAnExportProbe(),
	//buildSinusoidalInitialField(E,Z,{{1,1,0}},{{1.0,1.0,0.0}},Vector{{0.0,0.0,0.5}}),
	//SolverOptions{}
	//	.setTimeStep(5e-4)
	//	.setFinalTime(0.5)
	//	.setOrder(3)
	//};

	//GridFunction eOld{ solver.getFields().E[Z] };
	//auto normOld{ solver.getFields().getNorml2() };
	//solver.run();
	//GridFunction eNew{ solver.getFields().E[Z] };

	//EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	//EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

TEST_F(Solver3DTest, rotated_M45X_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(10,1,1,Element::HEXAHEDRON, 3.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongXAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 50;

	mfem::Vector center({ 1.5,0.5,0.5 });
	mfem::Vector polarization(3);
	rotateMinus45degAlongXAxis(unitVec(Z), polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.2, center, polarization, 1, Source::CartesianAngles({ M_PI_4,0.0,0.0 })),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, rotated_M45Y_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(10,1,1,Element::HEXAHEDRON, 3.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongYAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 50;

	mfem::Vector center({ 1.5,0.5,0.5 }), rotCenter(3);
	rotateMinus45degAlongYAxis(center, rotCenter);
	mfem::Vector polarization(3);
	rotateMinus45degAlongYAxis(unitVec(Z), polarization);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.2, rotCenter, polarization, 1, Source::CartesianAngles({ 0.0,M_PI_4,0.0 })),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, rotated_M45Z_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(10,1,1,Element::HEXAHEDRON, 3.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongZAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 50;

	mfem::Vector center({ 1.5,0.5,0.5 }), rotCenter(3);
	mfem::Vector polarization(3);
	rotateMinus45degAlongZAxis(unitVec(Z), polarization);
	rotateMinus45degAlongZAxis(center, rotCenter);

	maxwell::Solver solver{
	model,
	probes,
	buildGaussianInitialField(E, 0.2, rotCenter, polarization, 1, Source::CartesianAngles({ 0.0,0.0,M_PI_4 })),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, rotated_AllDir_centered_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(10,1,1,Element::HEXAHEDRON, 3.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongXAxis);
	mesh.Transform(rotateMinus45degAlongYAxis);
	mesh.Transform(rotateMinus45degAlongZAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 50;

	mfem::Vector center({ 1.5,0.5,0.5 }), rCenter1(3), rCenter2(3), rotCenter(3);
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
	buildGaussianInitialField(E, 0.2, rotCenter, polarization, 1, Source::CartesianAngles({ M_PI_4, M_PI_4, M_PI_4 })),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, rotated_AllDir_upwind_hexa_1dot5)
{

	Mesh mesh{ Mesh::MakeCartesian3D(10,1,1,Element::HEXAHEDRON, 3.0, 1.0, 1.0) };
	mesh.Transform(rotateMinus45degAlongXAxis);
	mesh.Transform(rotateMinus45degAlongYAxis);
	mesh.Transform(rotateMinus45degAlongZAxis);

	AttributeToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh,AttributeToMaterial{},attToBdr,AttributeToInteriorBoundary{} };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 50;

	mfem::Vector center({ 1.5,0.5,0.5 }), rCenter1(3), rCenter2(3), rotCenter(3);
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
	buildGaussianInitialField(E, 0.2, rotCenter, polarization, 1, Source::CartesianAngles({ M_PI_4, M_PI_4, M_PI_4 })),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(6.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, compare_SpectralToBase_centered) {

	Probes probes;

	maxwell::Solver solver{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
		.setOrder(3)
		.setSpectralEO()
	};

	maxwell::Solver solverSpectral{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(6.0)
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

TEST_F(Solver3DTest, compare_SpectralToBase_upwind) {

	Probes probes;

	maxwell::Solver solver{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(6.0)
		.setOrder(3)
		.setSpectralEO()
	};

	maxwell::Solver solverSpectral{
	buildModel(10,1,1, Element::Type::HEXAHEDRON, 3.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.2, mfem::Vector({1.5,0.5,0.5}), unitVec(Z)),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(6.0)
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
