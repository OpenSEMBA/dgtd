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

	AttributeToBoundary buildAttrToBdrMap3D(const BdrCond& bdr1, const BdrCond& bdr2, const BdrCond& bdr3, const BdrCond& bdr4, const BdrCond& bdr5, const BdrCond& bdr6)
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

	Probes buildProbesWithAnExportProbe()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}


};

TEST_F(Solver3DTest, 3D_centered_hexa_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(5,1,1, Element::Type::HEXAHEDRON, 5.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, Z, 0.7, mfem::Vector({2.5,0.5,0.5})),
	SolverOptions{}
		.setTimeStep(1e-3)
		.setCentered()
		.setFinalTime(5.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, 3D_upwind_hexa_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(5,1,1, Element::Type::HEXAHEDRON, 5.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, Z, 0.7, mfem::Vector({2.5,0.5,0.5})),
	SolverOptions{}
		.setTimeStep(1e-3)
		.setFinalTime(5.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver3DTest, 3D_centered_tetra_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(5,1,1, Element::Type::TETRAHEDRON, 5.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, Z, 0.7, mfem::Vector({2.5,0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(5.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver3DTest, 3D_upwind_tetra_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 50;

	maxwell::Solver solver{
	buildModel(5,1,1, Element::Type::TETRAHEDRON, 5.0, 1.0, 1.0, BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, Z, 0.7, mfem::Vector({2.5,0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(5.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}


TEST_F(Solver3DTest, squareBox_3D_centered_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5, 0.5}},
		PointProbe{H, Y, {0.0, 0.5, 0.5}}
	};
	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(7,1,1, Element::Type::HEXAHEDRON,1.0,1.0,1.0,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, Z, 0.1, mfem::Vector({0.5,0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}
TEST_F(Solver3DTest, DISABLED_box_pec_upwind_3D)
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
