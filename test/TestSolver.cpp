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

namespace HelperFunctions {
		
	mfem::Array<int> getH1LexOrder(const mfem::H1_FECollection* fec) {
		auto* fe = fec->FiniteElementForGeometry(mfem::Geometry::SEGMENT);
		const mfem::NodalFiniteElement* nodal_fe =
			dynamic_cast<const mfem::NodalFiniteElement*>(fe);
		mfem::Array<int> lexOrder = nodal_fe->GetLexicographicOrdering();
		return lexOrder;
	}
	mfem::SparseMatrix operatorToSparseMatrix(const mfem::Operator* op) {
		
		int width = op->Width();
		int height = op->Height();
		mfem::SparseMatrix res(height,height);

		mfem::Vector x(width), y(height);
		x = 0.0;

		for (int i = 0; i < width; i++)
		{
			x(i) = 1.0;
			op->Mult(x, y);
			for (int j = 0; j < height; j++)
			{
				if (y(j) != 0.0)
				{
					res.Add(i, j, y[j]);
				}
			}
			x(i) = 0.0;
		}

		res.Finalize();
		return res;
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
	int nx = 21; int ny = 21; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);
	//const char* mesh_file = "ref-square.mesh";
	//mfem::Mesh *readmesh = new mfem::Mesh(mesh_file, 1, 1);
	//mfem::Mesh mesh = mfem::Mesh(*readmesh);

	Solver::Options opts;
	opts.order = 2;
	opts.dt = 1e-3;
	opts.t_final = 0.05;
	opts.vis_steps = 1;
	AnalyticalFunctions::standingWaveModeX = 1;
	AnalyticalFunctions::standingWaveModeY = 1;
	
	Solver solver(opts, mesh);
	solver.getMesh().GetBoundingBox(
		AnalyticalFunctions::meshBoundingBoxMin, 
		AnalyticalFunctions::meshBoundingBoxMax);

	
	solver.setInitialElectricField(AnalyticalFunctions::gaussianFunction);
	//solver.run();
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

	TEST(DG, checkMassMatrix)
	{
		int order = 1;
		int dimension = 1;
		FiniteElementCollection* fec;
		FiniteElementSpace* fes;

		Mesh mesh = Mesh::MakeCartesian1D(1);
		fec = new H1_FECollection(order, dimension);
		fes = new FiniteElementSpace(&mesh, fec);
		
		BilinearForm massMatrix(fes);
		massMatrix.AddDomainIntegrator(new MassIntegrator);
		massMatrix.Assemble();
		massMatrix.Finalize();

		EXPECT_NEAR(2.0 / 6.0, massMatrix(0, 0), 1e-3);
		EXPECT_NEAR(1.0 / 6.0, massMatrix(0, 1), 1e-3);
		EXPECT_NEAR(1.0 / 6.0, massMatrix(1, 0), 1e-3);
		EXPECT_NEAR(2.0 / 6.0, massMatrix(1, 1), 1e-3);

	}

	TEST(DG, checkMassMatrixIsSameForH1andDG)
	{
		const int maxOrder = 5;
		const int dimension = 1;
		Mesh meshH1 = Mesh::MakeCartesian1D(1);
		Mesh meshDG = Mesh::MakeCartesian1D(1);

		for (int order = 3; order < maxOrder; order++) {

			std::cout << "Checking order: " << order << std::endl;

			auto fecH1 = new H1_FECollection(order, dimension, BasisType::ClosedUniform);
			auto fecDG = new DG_FECollection(order, dimension, BasisType::ClosedUniform);

			FiniteElementSpace* fesH1 = new FiniteElementSpace(
				&meshH1, fecH1);
			FiniteElementSpace* fesDG = new FiniteElementSpace(
				&meshDG, fecDG);

			auto lexOrderH1 = HelperFunctions::getH1LexOrder(fecH1);

			const Operator* rotatorH1 = fesH1->GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC);

			mfem::SparseMatrix rotatorMatrix = HelperFunctions::operatorToSparseMatrix(rotatorH1);
			
			rotatorMatrix.PrintMatlab(std::cout);

			//GridFunction H1Nodes(fesH1);
			//meshH1.GetNodes(H1Nodes);
			//GridFunction DGNodes(fesDG);
			//meshDG.GetNodes(DGNodes);

			//GridFunction newH1Nodes(fesH1);
			//rotatorH1->Mult(H1Nodes, newH1Nodes);

			//rotatorH1->PrintMatlab(std::cout);

			BilinearForm massMatrixH1(fesH1);
			massMatrixH1.AddDomainIntegrator(new MassIntegrator);
			massMatrixH1.Assemble();
			massMatrixH1.Finalize();
			BilinearForm massMatrixDG(fesDG);
			massMatrixDG.AddDomainIntegrator(new MassIntegrator);
			massMatrixDG.Assemble();
			massMatrixDG.Finalize();
			BilinearForm rotatedMassMatrixH1(fesH1);
			rotatedMassMatrixH1.AddDomainIntegrator(new MassIntegrator);
			rotatedMassMatrixH1.Assemble();
			rotatedMassMatrixH1.Finalize();


			//massMatrixH1.Print(std::cout);
			//std::cout << std::endl;
			//massMatrixDG.Print(std::cout);
			//std::cout << std::endl;

			ASSERT_EQ(massMatrixH1.NumRows(), massMatrixDG.NumRows());
			ASSERT_EQ(massMatrixH1.NumCols(), massMatrixDG.NumCols());

			for (std::size_t i = 0; i < massMatrixDG.NumRows(); i++) {
				for (std::size_t j = 0; j < massMatrixDG.NumCols(); j++) {
					EXPECT_NEAR(massMatrixDG(i, j), massMatrixH1(i, j), 1e-5);
				}
			}
		}
	}
}
