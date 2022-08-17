#include "gtest/gtest.h"

#include "AnalyticalFunctions2D.h"
#include "maxwell/Solver.h"

using namespace maxwell;
using namespace mfem;
using namespace AnalyticalFunctions2D;

class TestSolver2D : public ::testing::Test {
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

TEST_F(TestSolver2D, DISABLED_checkRun)
{
	/*The purpose of this test is to check the run() function for the solver object
	and test the different available options.
	
	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.
	
	Then, a solver object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction ez_)
	
	Lastly, the run() function is called.*/

	
	int nx = 101; int ny = 1; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges,100,1);

	maxwell::Solver::Options opts;
	opts.order = 2;
	opts.dt = 5e-4;
	opts.t_final = 1000*opts.dt;
	opts.vis_steps = 500;

	standingWaveModeX = 1;
	standingWaveModeY = 1;
	
	Solver solver(opts, mesh);
	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);

	
	solver.setInitialElectricField(gaussianFunction);
	solver.run();
}
TEST_F(TestSolver2D, DISABLED_checkMeshDimensions) 
{
	/*This test ensures that the number of elements of any 2D Cartesian
	mesh is equal to the product of the horizontal and vertical segments
	
	Dimensional parameters are declared which are then used to create
	a mesh object, then the test comparison is made.*/


	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	EXPECT_EQ(nx*ny, mesh.GetNE());

}
TEST_F(TestSolver2D, DISABLED_checkMeshElementVertices) 
{
	/*This test was created to understand the process of mesh creation
	and assignation of vertex index to elements.
	
	First, dimensional variables are declared, which are then used
	to create a new mesh object.
	
	Then firstElementVerticesVector and lastElementVerticesVector are 
	initialized and assigned values manually, with the values we expect
	these elements will have. We also retrieve the vertices for the first
	and last element of the mesh, and then store them inside vectors.
	
	Lastly, we compare that the vertices we retrieved from the mesh are
	equal to those we presumed at the start.*/

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	std::vector<int> firstElementVerticesVector = { 0, 1, nx+2, nx+1 };
	std::vector<int> lastElementVerticesVector = { nx-1, nx, nx*2+1, nx*2 };
	mfem::Array<int> meshArrayFirstElement;
	mfem::Array<int> meshArrayLastElement;

	mesh.GetElementVertices(0, meshArrayFirstElement);
	mesh.GetElementVertices(nx*ny-1, meshArrayLastElement);

	std::vector<int> vectorFirstElement(meshArrayFirstElement.begin(), meshArrayFirstElement.end());
	std::vector<int> vectorLastElement(meshArrayLastElement.begin(), meshArrayLastElement.end());

	EXPECT_EQ(firstElementVerticesVector, vectorFirstElement);
	EXPECT_EQ(lastElementVerticesVector, vectorLastElement);

}
TEST_F(TestSolver2D, DISABLED_checkMeshBoundaries)
{
	/*In this test we aim to compare the boundary DoFs for a mesh generated through
	the Mesh class constructors for a 2DCartesian and 'square3x3.mesh', a handcrafted mesh.
	
	First, we declare and initialise the variables we will require in order to create a
	FiniteElementSpace object, that is a FiniteElementCollection and a Mesh. In this we create
	two FES, one for the Cartesian2D mesh (Auto) and one for the handcrafted mesh (Manual).

	We then extract the Boundary DoFs from each of the FES.

	Lastly, we compare both of the Arrays where the DoFs are stored, to confirm that when it comes
	to boundaries of this nature, the Cartesian2D constructor works equally to the handcrafted mesh.
	*/

	int order = 1;
	int dimension = 2;
	int nx = 3; int ny = 3; bool generateEdges = true;
	mfem::Mesh meshAuto = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	const char* mesh_file = "square3x3.mesh";
	mfem::Mesh* readmesh = new mfem::Mesh(mesh_file, 1, 1);
	mfem::Mesh meshManual = mfem::Mesh(*readmesh);

	auto fec = new mfem::DG_FECollection(order, dimension, mfem::BasisType::GaussLegendre);

	auto fesAuto = new mfem::FiniteElementSpace(&meshAuto, fec);
	auto fesManual = new mfem::FiniteElementSpace(&meshManual, fec);

	mfem::Array<int> ess_tdof_list_auto;
	if (meshAuto.bdr_attributes.Size())
	{
		mfem::Array<int> ess_bdr_auto(meshAuto.bdr_attributes.Max());
		ess_bdr_auto = 1;
		fesAuto->GetEssentialTrueDofs(ess_bdr_auto, ess_tdof_list_auto);
	}

	mfem::Array<int> ess_tdof_list_manual;
	if (meshManual.bdr_attributes.Size())
	{
		mfem::Array<int> ess_bdr_manual(meshManual.bdr_attributes.Max());
		ess_bdr_manual = 1;
		fesManual->GetEssentialTrueDofs(ess_bdr_manual, ess_tdof_list_manual);
	}

	EXPECT_EQ(ess_tdof_list_auto, ess_tdof_list_manual);
}
TEST_F(TestSolver2D, DISABLED_mapMeshElementAndVertex) 
{

	/* This test was created with the aim to understand the mapping and ordering process
	of a mesh in a more visual way. It uses the mapQuadElementTopLeftVertex() function
	which, for a Quadrilateral Element, it extracts its top left vertex, which allows for a nigh
	full mapping of the mesh.
	
	First, dimensional variables are declared and a mesh is constructed.
	
	Then, the mapQuadElementTopLeftVertex extracts the top left vertex of each element and stores them
	in an integer vector.
	
	Lastly, we compare that the first mapped vertex is the first created vertex in the mesh 0,
	the top left vertex for the uppermost, rightmost element is equal to the last element's index - 1 
	(due to how mesh mapping works), and the size of the mapped vertices vector is equal to the number 
	of elements in the mesh - 1 (as it starts with index 0).*/

	int nx = 5; int ny = 5; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	std::vector<int> mapped = mapQuadElementTopLeftVertex(mesh);

	EXPECT_EQ(0, mapped[0]);
	EXPECT_EQ(nx - 1, mapped[mapped.size() -1]);
	EXPECT_EQ(nx*ny-1, mapped.size()-1);
	
}
TEST_F(TestSolver2D, DISABLED_centered_flux_AMR)
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

	const char* mesh_file = "amr-quad.mesh";
	Mesh mesh(mesh_file);
	mesh.UniformRefinement();
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Sources sources;
	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 })));

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(2.92);
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(model, Probes(), sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);

	EXPECT_GT(eOld.Max(), eNew.Max());
}
TEST_F(TestSolver2D, DISABLED_checkMeshInvariance) 
{
	/* This test aimed to prove mesh invariance after being used by our Solver object, and
	* deep the understanding of C++ programming.
	
	First, dimensional variables are declared, a mesh object is constructed and copied for comparison purposes,
	a Solver object is created by inputting the mesh as a parameter.
	
	Then, the copied mesh and solver object's mesh get their top left vertex extracted for each element.
	
	Lastly, several comparisons are made to the copied mesh and the solver's mesh, such as dimension, vertex index
	and size.*/

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	auto meshCopy = mfem::Mesh(mesh);
	Solver2D solver(Solver2D::Options(), mesh);

	std::vector<int> meshMap = mapQuadElementTopLeftVertex(meshCopy);
	std::vector<int> solverMeshMap = mapQuadElementTopLeftVertex(solver.getMesh());

	ASSERT_EQ(meshCopy.Dimension(), solver.getMesh().Dimension());
	EXPECT_EQ(meshMap[0], solverMeshMap[0]);
	EXPECT_EQ(meshMap.size(), solverMeshMap.size());
}
TEST_F(TestSolver2D, DISABLED_periodic) //TODO ADD ENERGY CHECK
{
	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
	Vector periodic({ 0.0, 1.0 });
	std::vector<Vector> trans;
	trans.push_back(periodic);
	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D, mesh2D.CreatePeriodicVertexMapping(trans));

	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());

	Probes probes;
	//probes.addExporterProbeToCollection(ExporterProbe());
	//probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.2, 0.0 })));

	maxwell::Solver solver(model, probes, sources, buildDefaultSolverOpts(1.0));

	solver.run();

}
TEST_F(TestSolver2D, DISABLED_periodic_strong) //TODO ADD ENERGY CHECK
{
	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
	Vector periodic({ 0.0, 1.0 });
	std::vector<Vector> trans;
	trans.push_back(periodic);
	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D, mesh2D.CreatePeriodicVertexMapping(trans));

	maxwell::Solver::Options opts;
	opts.evolutionOperatorOptions = FiniteElementEvolution::Options();
	opts.evolutionOperatorOptions.disForm = DisForm::Strong;

	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());

	Probes probes;
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.0, 0.0 })));

	maxwell::Solver solver(
		model,
		probes,
		sources,
		opts);

	solver.run();

}
TEST_F(TestSolver2D, DISABLED_centered_NC_MESH) //TODO ADD ENERGY CHECK
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

	const char* mesh_file = "star-mixed.mesh";
	Mesh mesh(mesh_file);
	mesh.UniformRefinement();
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	Probes probes;
	//probes.addExporterProbeToCollection(ExporterProbe());
	//probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 })));

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(2.92);
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(model, probes, sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);

	EXPECT_GT(eOld.Max(), eNew.Max());
}
TEST_F(TestSolver2D, DISABLED_twoDimensionalResonantBox)
{
	Mesh mesh2D = Mesh::MakeCartesian2D(21, 21, Element::Type::QUADRILATERAL);
	std::vector<Attribute> attArrSingle = std::vector<Attribute>({ 1 });
	Material mat11 = Material(1.0, 1.0);
	std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
	AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatMap(attArrSingle, matArrSimple);
	AttributeToBoundary attToBdrVec;
	Model model(mesh2D, attToMatVec, attToBdrVec);

	double spread = 2.0;
	double coeff = 20.0;
	const Vector dev = Vector({ 0.0,0.0 });
	Source EXFieldSource = Source(model, spread, coeff, dev, X, E); 
	Sources sources;
	sources.addSourceToVector(EXFieldSource);

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 100;

	maxwell::Solver::Options solverOpts;

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.dt = 1e-4;
	solverOpts.order = 1;

	maxwell::Solver solver(model, probes,
		sources, solverOpts);

	solver.run();

}