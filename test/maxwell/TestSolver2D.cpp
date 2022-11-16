#include "gtest/gtest.h"

#include "AnalyticalFunctions2D.h"
#include "SourceFixtures.h"
#include "maxwell/Solver.h"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;
using namespace AnalyticalFunctions2D;

class TestSolver2D : public ::testing::Test {
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

	AttributeToBoundary buildAttrToBdrMap2D(const BdrCond& bdrB, const BdrCond& bdrR, const BdrCond& bdrT, const BdrCond& bdrL)
	{
		return {
			{1, bdrB},
			{2, bdrR},
			{3, bdrT},
			{4, bdrL},
		};
	}

	Probes buildExportProbes()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}

};

TEST_F(TestSolver2D, box_pec_centered_2D)
{
	/*The purpose of this test is to check the run() function for the solver object
	and test the different available options.
	
	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.
	
	Then, a solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction ez_)
	
	Lastly, the run() function is called.*/

	maxwell::Solver solver{
	buildModel(5,5),
	buildExportProbes(),
	buildGaussianInitialField(E, Z, 0.1, 0.5, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(0.5)
		.setOrder(3)
	};

	GridFunction eOld{ solver.getFields().E[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Z] };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

TEST_F(TestSolver2D, box_pec_centered_square_2D)
{
	/*The purpose of this test is to check the run() function for the solver object
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction ez_)

	Lastly, the run() function is called.*/

	maxwell::Solver solver{
	buildModel(5,5,Element::Type::QUADRILATERAL),
	buildExportProbes(),
	buildGaussianInitialField(E, Z, 0.1, 0.5, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setCentered()
		.setFinalTime(0.5)
		.setOrder(4)
	};

	GridFunction eOld{ solver.getFields().E[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Z] };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

TEST_F(TestSolver2D, box_pec_upwind_2D)
{
	/*The purpose of this test is to check the run() function for the solver object
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction ez_)

	Lastly, the run() function is called.*/

	maxwell::Solver solver{
	buildModel(5,5),
	buildExportProbes(),
	buildGaussianInitialField(E, Z, 0.1, 0.5, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(0.5)
		.setOrder(3)
	};

	GridFunction eOld{ solver.getFields().E[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Z] };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

TEST_F(TestSolver2D, box_pec_upwind_square_2D)
{
	/*The purpose of this test is to check the run() function for the solver object
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction ez_)

	Lastly, the run() function is called.*/

	maxwell::Solver solver{
	buildModel(5,5,Element::Type::QUADRILATERAL),
	buildExportProbes(),
	buildGaussianInitialField(E, Z, 0.1, 0.5, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(0.5)
		.setOrder(4)
	};

	GridFunction eOld{ solver.getFields().E[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Z] };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

//TEST_F(TestSolver2D, DISABLED_centered_flux_AMR)
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
//TEST_F(TestSolver2D, DISABLED_periodic) //TODO ADD ENERGY CHECK
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
//TEST_F(TestSolver2D, DISABLED_periodic_strong) //TODO ADD ENERGY CHECK
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
//TEST_F(TestSolver2D, DISABLED_centered_NC_MESH) //TODO ADD ENERGY CHECK
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
//TEST_F(TestSolver2D, DISABLED_resonantBox)
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