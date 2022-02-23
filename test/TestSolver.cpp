#include "gtest/gtest.h"
#include <math.h>

#include "Solver.h"

using namespace Maxwell;

namespace AnalyticalFunctions {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;
	std::size_t standingWaveModeX = 1, standingWaveModeY = 1;

	const double PI = atan(1.0) * 4;

	double gaussianFunction(const Solver::Position& pos)
	{
		mfem::Vector normalizedPos(2);
		for (size_t i = 0; i < 2; i++) {
			double center = (meshBoundingBoxMin[i] + meshBoundingBoxMax[i]) * 0.5;
			normalizedPos[i] = 2 * (pos[i] - center) / (meshBoundingBoxMax[i] - meshBoundingBoxMin[i]);
		}

		return exp(-20. * (pow(normalizedPos[0], 2) + pow(normalizedPos[1], 2)));
	}

	double gaussianFunction1D(const mfem::Vector pos)
	{
		mfem::Vector normalizedPos(2);
		for (size_t i = 0; i < 2; i++) {
			double center = 0.5;
			normalizedPos[i] = 2 * (pos[i] - center);
		}

		return exp(-20. * (pow(normalizedPos[0], 2)));
	}

	double standingWaveFunction(const Solver::Position& pos)
	{
		mfem::Vector normalizedPos(2);
		mfem::Vector L(2);
		for (size_t i = 0; i < 2; i++) {
			double center = (meshBoundingBoxMin[i] + meshBoundingBoxMax[i]) * 0.5;
			L[i] = meshBoundingBoxMax[i] - meshBoundingBoxMin[i];
			normalizedPos[i] = (pos[i] - meshBoundingBoxMin[i]) / L[i];
		}

		return sin(normalizedPos[0]  * PI * double(standingWaveModeX)) *
			sin(normalizedPos[1] * PI * double(standingWaveModeY));
	}
}

class TestSolver : public ::testing::Test {
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
TEST_F(TestSolver, checkRun)
{
	int nx = 9; int ny = 9; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	//const char* mesh_file = "square3x3.mesh";
	//mfem::Mesh *readmesh = new mfem::Mesh(mesh_file, 1, 1);
	//mfem::Mesh mesh = mfem::Mesh(*readmesh);

	Solver::Options opts;
	opts.order = 2;
	//opts.dt = 1e-3;
	opts.t_final = 1;
	opts.vis_steps = 1;
	AnalyticalFunctions::standingWaveModeX = 1;
	AnalyticalFunctions::standingWaveModeY = 1;
	
	Solver solver(opts, mesh);
	solver.getMesh().GetBoundingBox(
		AnalyticalFunctions::meshBoundingBoxMin, 
		AnalyticalFunctions::meshBoundingBoxMax);

	
	solver.setInitialElectricField(AnalyticalFunctions::gaussianFunction);
	solver.run();
}
TEST_F(TestSolver, checkMeshDimensions) 
{

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	EXPECT_EQ(nx*ny, mesh.GetNE());

}
TEST_F(TestSolver, checkMeshElementVertices) 
{

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
TEST_F(TestSolver, checkMeshBoundaries)
{
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

	meshAuto.Print(std::cout);
	std::cout << std::endl;
	meshManual.Print(std::cout);
	std::cout << std::endl;

	EXPECT_EQ(ess_tdof_list_auto, ess_tdof_list_manual);
}
TEST_F(TestSolver, mapMeshElementAndVertex) 
{

	int nx = 5; int ny = 5; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	std::vector<int> mapped = mapQuadElementTopLeftVertex(mesh);

	EXPECT_EQ(0, mapped[0]);
	EXPECT_EQ(nx*ny-1, mapped.size()-1);
	EXPECT_EQ(nx-1, mapped[mapped.size()-1]);
}
TEST_F(TestSolver, checkMeshInvariance) 
{

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	Solver solver(Solver::Options(), mesh);

	std::vector<int> meshMap = mapQuadElementTopLeftVertex(mesh);
	std::vector<int> solverMeshMap = mapQuadElementTopLeftVertex(solver.getMesh());

	ASSERT_EQ(mesh.Dimension(), solver.getMesh().Dimension());
	EXPECT_EQ(meshMap[0], solverMeshMap[0]);
	EXPECT_EQ(meshMap.size(), solverMeshMap.size());
}

namespace mfem {
TEST_F(TestSolver, oneDimensional)
{
	int order = 1;
	const int dimension = 1;
	FiniteElementCollection* fec;
	FiniteElementSpace* fes;

	Mesh mesh = Mesh::MakeCartesian1D(1);
	fec = new DG_FECollection(order, dimension, BasisType::GaussLobatto);
	fes = new FiniteElementSpace(&mesh, fec);

	BilinearForm massMatrix(fes);
	massMatrix.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	massMatrix.Assemble();
	massMatrix.Finalize();

	ConstantCoefficient one(1.0);
	BilinearForm K(fes);
	K.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
	K.Assemble();
	K.Finalize();

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	BilinearForm fluxForm(fes);
	fluxForm.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
	fluxForm.AddInteriorFaceIntegrator(new DGTraceIntegrator(n[0], 0.0, 1.0));
	fluxForm.Assemble();
	fluxForm.Finalize();

	double gaussianFunction = exp(-20. * (pow(0, 2) + pow(1, 2)));

	GridFunction Ez(fes);
	Ez.ProjectCoefficient(FunctionCoefficient(AnalyticalFunctions::gaussianFunction1D));
	GridFunction Hy(fes);

	std::unique_ptr<ParaViewDataCollection> pd = NULL;
	pd = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh);
	pd->SetPrefixPath("ParaView");
	pd->RegisterField("ez", &Ez);
	pd->RegisterField("hy", &Hy);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);

	double time = 0.0;
	bool done = false;
	double dt = 1e-2;
	double t_final = 2.0;
	int vis_steps = 1;

	Vector aux(fes->GetVSize());
	Vector ezNew(fes->GetVSize());
	Vector hyNew(fes->GetVSize());

	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();

	for (int cycle = 0; !done;)
	{

		// Update E.
		K.Mult(Hy, aux);
		massMatrix.Mult(aux, ezNew);
		ezNew *= -dt;
		ezNew.Add(1.0, Ez);

		// Update H.
		K.Mult(ezNew, aux);
		massMatrix.Mult(aux, hyNew);
		hyNew *= -dt;
		hyNew.Add(1.0, Hy);

		Ez = ezNew;
		Hy = hyNew;

		time += dt;
		cycle++;

		done = (time >= t_final - 1e-8 * dt);

		if (done || cycle % vis_steps == 0) {
			pd->SetCycle(cycle);
			pd->SetTime(time);
			pd->Save();
		}
	}
	

}

}

