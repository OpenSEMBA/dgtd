#include "gtest/gtest.h"

#include "AnalyticalFunctions2D.h"
#include "SourceFixtures.h"
#include "maxwell/Solver.h"

#include "maxwell/Types.h"
#include "mfem.hpp"
#include "maxwell/Model.h"
#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/mfemExtension/BilinearForm_IBFI.hpp"
#include "maxwell/mfemExtension/LinearIntegrators.h"
#include "maxwell/mfemExtension/LinearForm_IBFI.hpp"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;
using namespace AnalyticalFunctions2D;

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

		return Model(Mesh::MakeCartesian2D(nx,ny,elType), AttributeToMaterial{}, buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL));
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

		return Model(Mesh::MakeCartesian2D(nx, ny, elType, false, sx, sy), AttributeToMaterial{}, buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL));
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

	Probes buildProbesWithAnExportProbe()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}

	Vector fieldCenter{ { 0.5, 0.5 } };
	Source::Polarization zPolarization()
	{
		return Source::Polarization({ 0.0, 0.0, 1.0 });
	}

};

TEST_F(Solver2DTest, 2D_pec_centered_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	}; 

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//At the left boundary the electric field should be closed to zero and
	//the magnetic field reaches a maximum close to 1.0 
	//(the wave splits in two and doubles at the boundary).
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver2DTest, 2D_pec_centered_quadrilaterals_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(), 1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 
	// (the wave splits in two and doubles at the boundary).
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver2DTest, 2D_pec_upwind_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 
	// (the wave splits in two and doubles at the boundary).
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver2DTest, 2D_pec_upwind_quadrilaterals_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 
	// (the wave splits in two and doubles at the boundary).
	auto eMaxFrameMid{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrameMid.second, tolerance);

	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

}

TEST_F(Solver2DTest, 2D_sma_upwind_quadrilaterals_1dot5D)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(10, 1, mfem::Element::Type::QUADRILATERAL,1.0, 1.0, BdrCond::PMC, BdrCond::SMA, BdrCond::PMC, BdrCond::SMA),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	GridFunction eOld{ solver.getFields().E[Z] };

	auto zeros{ eOld };
	zeros = 0.0;
	EXPECT_TRUE( eOld.DistanceTo(zeros) > 1e-2);

	solver.run();

	double tolerance{ 1e-2 };
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMin() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMin() };
	EXPECT_NEAR(0.0, hMaxFrame.second, tolerance);


}

TEST_F(Solver2DTest, 2D_rotated_quadrilateral_centered_1dot5D)
{
	auto mesh{ Mesh::LoadFromFile("./testData/severalrotatedquads.mesh",1,0) };
	mesh.UniformRefinement();
	AttributeToBoundary attToBdr{ {1,BdrCond::PEC}, {2,BdrCond::PMC}};
	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorBoundary{} };
	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 30;

	fieldCenter = Vector({ 2.0, 2.0 });

	maxwell::Solver solver {
		model,
		probes,
		buildGaussianInitialField(E, 0.5, fieldCenter, zPolarization(), 1, Source::CartesianAngles({0.0,0.0,-M_PI_4})),
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

TEST_F(Solver2DTest, 2D_pec_centered_1dot5D_spectral)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
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
	// the magnetic field reaches a maximum close to 1.0 
	// (the wave splits in two and doubles at the boundary).
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver2DTest, 2D_pec_upwind_1dot5D_spectral)
{

	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}}
	};

	probes.exporterProbes[0].visSteps = 40;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 
	// (the wave splits in two and doubles at the boundary).
	auto eMaxFrame{ solver.getPointProbe(0).findFrameWithMax() };
	EXPECT_NEAR(0.0, eMaxFrame.second, tolerance);
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver2DTest, compare_2DSpectralToBase_centered) {

	Probes probes;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	};

	maxwell::Solver solverSpectral{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	Vector zeroVec{ solver.getFields().E[Z].Size() };
	zeroVec = 0.0;
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().E[Z].Elem(i) - solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[Y].Elem(i) - solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[X].Elem(i) - solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().E[Z].Elem(i) - solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[Y].Elem(i) - solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[X].Elem(i) - solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}

}

TEST_F(Solver2DTest, compare_2DSpectralToBase_upwind) {

	Probes probes;

	maxwell::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	maxwell::Solver solverSpectral{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, fieldCenter, zPolarization(),1),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(3)
		.setSpectralEO()
	};

	Vector zeroVec{ solver.getFields().E[Z].Size() };
	zeroVec = 0.0;
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().E[Z].Elem(i) - solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[Y].Elem(i) - solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[X].Elem(i) - solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < zeroVec.Size(); ++i) {
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().E[Z].Elem(i) - solverSpectral.getFields().E[Z].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[Y].Elem(i) - solverSpectral.getFields().H[Y].Elem(i), 1e-5);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getFields().H[X].Elem(i) - solverSpectral.getFields().H[X].Elem(i), 1e-5);
	}

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

//TEST_F(Solver2DTest, DISABLED_centered_flux_AMR)
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
//	const char* mesh_file = "amr-quad.mesh";
//	Mesh mesh(mesh_file);
//	mesh.UniformRefinement();
//	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());
//
//	Sources sources;
//	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 })));
//
//	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(2.92);
//	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
//
//	maxwell::Solver solver(model, Probes(), sources, solverOpts);
//
//	GridFunction eOld = solver.getFieldInDirection(E, Z);
//	solver.run();
//	GridFunction eNew = solver.getFieldInDirection(E, Z);
//
//	EXPECT_GT(eOld.Max(), eNew.Max());
//}
//TEST_F(Solver2DTest, DISABLED_periodic) //TODO ADD ENERGY CHECK
//{
//	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
//	Vector periodic({ 0.0, 1.0 });
//	std::vector<Vector> trans;
//	trans.push_back(periodic);
//	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D, mesh2D.CreatePeriodicVertexMapping(trans));
//
//	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());
//
//	Probes probes;
//	//probes.addExporterProbeToCollection(ExporterProbe());
//	//probes.vis_steps = 20;
//
//	Sources sources;
//	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.2, 0.0 })));
//
//	maxwell::Solver solver(model, probes, sources, buildDefaultSolverOpts(1.0));
//
//	solver.run();
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
//	opts.evolutionOperatorOptions = FiniteElementEvolution::Options();
//	opts.evolutionOperatorOptions.disForm = DisForm::Strong;
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
//	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
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
//	solverOpts.dt = 1e-4;
//	solverOpts.order = 1;
//
//	maxwell::Solver solver(model, probes,
//		sources, solverOpts);
//
//	solver.run();
//
//}