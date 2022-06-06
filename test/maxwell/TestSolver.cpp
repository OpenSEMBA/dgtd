#include "gtest/gtest.h"
#include <math.h>

#include "maxwell/Solver.h"

using namespace maxwell;

using Interval = std::pair<double, double>;

namespace AnalyticalFunctions1D {
	mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;

	const double PI = atan(1.0) * 4;

	double gaussianFunction(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 2.0 * (pos[0] - center) /
			            ((meshBoundingBoxMax[0] - meshBoundingBoxMin[0]));
		
		return exp(-20. * pow(normalizedPos, 2));
	}

	double gaussianFunctionHalfWidth(const mfem::Vector pos)
	{
		double normalizedPos;
		double center = (meshBoundingBoxMin[0] + meshBoundingBoxMax[0]) * 0.5;
		normalizedPos = 4.0 * (pos[0] - center/2) /
			((meshBoundingBoxMax[0] - meshBoundingBoxMin[0]));

		return exp(-20. * pow(normalizedPos, 2));
	}
}

namespace HelperFunctions {

	void setAttributeIntervalMesh1D(
		const std::map<Attribute,Interval>& attToInterval,
		Mesh& mesh)
	{
		for (auto const& kv : attToInterval) {
			DenseMatrix changeAttMat(1, 2);
			changeAttMat.Elem(0, 0) = kv.second.first;
			changeAttMat.Elem(0, 1) = kv.second.second;
			Array<int> elemID;
			Array<IntegrationPoint> integPoint;
			mesh.FindPoints(changeAttMat, elemID, integPoint);

			if (elemID.begin() > elemID.end()) {
				throw std::exception("Lower Index bigger than Higher Index.");
			}
			if (elemID[1] > mesh.GetNE()) {
				throw std::exception("Declared element index bigger than Mesh Number of Elements.");
			}
			for (int i = elemID[0]; i <= elemID[1]; i++) {
				mesh.SetAttribute(i, kv.first);
			}
		}
	}

	std::vector<int> mapQuadElementTopLeftVertex(
		const mfem::Mesh& mesh)
	{
		std::vector<int> res;
		for (int i = 0; i < mesh.GetNE(); i++) {
			mfem::Array<int> meshArrayElement;
			mesh.GetElementVertices(i, meshArrayElement);
			res.push_back(meshArrayElement[0]);
		}

		return res;
	}

	std::map<Time, FieldFrame>::const_iterator findTimeId(
		const std::map<Time, FieldFrame>& timeMap,
		const Time& timeToFind,
		const double tolerance)
	{
		for (auto it = timeMap.begin(); it != timeMap.end(); it++) {
			const Time& time = it->first;
			if (abs(time - timeToFind) < tolerance) {
				return it;
			}
		}
		return timeMap.end();
	}

}

using namespace AnalyticalFunctions1D;
using namespace HelperFunctions;

class TestMaxwellSolver : public ::testing::Test {
protected:

	Model buildOneDimOneMatModel(
		const int meshIntervals = 51, 
		const BdrCond& bdrL = BdrCond::PEC, 
		const BdrCond& bdrR = BdrCond::PEC) {

		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), AttributeToMaterial(), buildAttrToBdrMap1D(bdrL,bdrR));
	}

	Sources buildSourcesWithDefaultSource(
		Model& model, 
		const FieldType& ft = E,
		const Direction& d = X, 
		const double spread = 2.0, 
		const double coeff = 1.0, 
		const Vector dev = Vector({ 0.0 })) {

		Sources res;
		res.addSourceToVector(Source(model, ft, d,  spread, coeff, dev));
		return res;
	}

	Probes buildProbesWithDefaultPointsProbe(
		const FieldType& fToExtract = E, 
		const Direction& dirToExtract = X)
	{
		Probes res;
		res.vis_steps = 20;
		res.extractDataAtPoints = true; 
		res.addProbeToCollection(PointsProbe(fToExtract, dirToExtract,
			std::vector<std::vector<double>>{{0.0},{0.5},{1.0}}));
		return res;
	}

	maxwell::Solver::Options buildDefaultSolverOpts(const double tFinal = 2.0)
	{
		maxwell::Solver::Options res;

		res.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
		res.t_final = tFinal;

		return res;
	}

	AttributeToBoundary buildAttrToBdrMap1D(const BdrCond& bdrL, const BdrCond& bdrR)
	{
		return {
			{1, bdrL},
			{2, bdrR}
		};
	}

	AttributeToMaterial buildAttToVaccumOneMatMap1D()
	{
		return { 
			{ 1, Material(1.0, 1.0) } 
		};
	}

	double getBoundaryFieldValueAtTime(
		const PointsProbe& probe,
		const Time& timeToFind,
		const int denseMatPointByOrder,
		const Direction& d)
	{
		auto itpos = HelperFunctions::findTimeId(probe.getConstFieldMovie(), timeToFind, 1e-6);
		if (itpos == probe.getConstFieldMovie().end()) {
			throw std::exception("Time value has not been found within the specified tolerance.");
		}
		auto FieldValueForTimeAtPoint = itpos->second.at(denseMatPointByOrder).at(d);

		return FieldValueForTimeAtPoint;
	}

};
TEST_F(TestMaxwellSolver, oneDimensional_centered)
{	
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts, with its mesh being one-dimensional.
	The field along Ey is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that after two seconds with PEC boundary conditions, the wave evolves
	back to its initial state within the specified error.*/

	Model model = buildOneDimOneMatModel();

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(
		model, 
		Probes(),
		buildSourcesWithDefaultSource(model), 
		solverOpts);
	
	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_centered_energy)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts, with its mesh being one-dimensional.
	The field along Ey is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that after two seconds with PEC boundary conditions, the wave evolves
	back to its initial state within the specified error.*/

	Model model = buildOneDimOneMatModel();

	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;

	maxwell::Solver solver(
		model,
		Probes(),
		buildSourcesWithDefaultSource(model),
		solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	EXPECT_GE(pow(eOld.Norml2(),2.0) + pow(hOld.Norml2(),2.0), pow(eNew.Norml2(),2.0) + pow(hNew.Norml2(),2.0));

}

TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC_EX)
{

	Model model = buildOneDimOneMatModel();

	maxwell::Solver solver(
		model, 
		buildProbesWithDefaultPointsProbe(),
		buildSourcesWithDefaultSource(model),
		buildDefaultSolverOpts());

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC_EY)
{
	Model model = buildOneDimOneMatModel();

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(E, Y),
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts());

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0, Y), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2, Y), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0, Y), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2, Y), 2e-3);

	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1, Y));
	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1, Y));

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC_EZ)
{

	Model model = buildOneDimOneMatModel();

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(E, Z),
		buildSourcesWithDefaultSource(model, E, Z),
		buildDefaultSolverOpts());

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0, Z), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2, Z), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0, Z), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2, Z), 2e-3);

	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1, Z));
	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1, Z));

}


TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC_HX)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::PMC, BdrCond::PMC);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(H, X),
		buildSourcesWithDefaultSource(model, H, X),
		buildDefaultSolverOpts());

	GridFunction hOld = solver.getFieldInDirection(H, X);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, X);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC_HY)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::PMC, BdrCond::PMC);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(H, Y),
		buildSourcesWithDefaultSource(model, H, Y),
		buildDefaultSolverOpts()
	);

	GridFunction hOld = solver.getFieldInDirection(H, Y);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Y);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0, Y), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2, Y), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0, Y), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2, Y), 2e-3);

	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1, Y));
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1, Y));

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC_HZ)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::PMC, BdrCond::PMC);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(H, Z),
		buildSourcesWithDefaultSource(model, H, Z),
		buildDefaultSolverOpts());

	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0, Z), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2, Z), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0, Z), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2, Z), 2e-3);

	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1, Z));
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1, Z));

}

TEST_F(TestMaxwellSolver, oneDimensional_upwind_SMA_EX)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(E, X),
		buildSourcesWithDefaultSource(model, E, X),
		buildDefaultSolverOpts(1.0));

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_SMA_EY)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	auto probes = buildProbesWithDefaultPointsProbe(E, Y);
	auto probeZ = PointsProbe(E, Z, std::vector<std::vector<double>>({ {0.0},{0.5},{1.0} }));
	probes.addProbeToCollection(probeZ);
	probes.addProbeToCollection(ExporterProbe());
	

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts(1.0));

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 0, Y), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 1, Y), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 2, Y), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 0, Z), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 1, Z), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 2, Z), 2e-3);

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_SMA_EZ)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(E, Z),
		buildSourcesWithDefaultSource(model, E, Z),
		buildDefaultSolverOpts(1.0));

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_GE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0, Z));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1, Z));
	EXPECT_GE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2, Z));

}

TEST_F(TestMaxwellSolver, twoSourceWaveTravelsToTheRight_SMA)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	Probes probes;
	probes.extractDataAtPoints = true;
	probes.addProbeToCollection(PointsProbe(E, Y, std::vector<std::vector<double>>{ {0.5}, { 0.8 } }));

	Sources sources;
	sources.addSourceToVector(Source(model, E, Y, 2.0, 1.0, Vector({ 0.0 })));
	sources.addSourceToVector(Source(model, H, Z, 2.0, 1.0, Vector({ 0.0 })));

	maxwell::Solver solver(
		model,
		probes,
		sources,
		buildDefaultSolverOpts(0.7));

	solver.run();

	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.3, 1, Y),
				getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0, Y),
				2e-3);

}
TEST_F(TestMaxwellSolver, twoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	Mesh mesh1D = Mesh::MakeCartesian1D(101);
	setAttributeIntervalMesh1D({ { 2, std::make_pair(0.76, 1.0) } }, mesh1D);

	Model model = Model(
		mesh1D, 
		{
			{1, Material(1.0, 1.0)},
			{2, Material(2.0, 1.0)}
		},
		{
			{1, BdrCond::SMA},
			{2, BdrCond::PEC}
		}
	);

	Probes probes;
	probes.extractDataAtPoints = true;
	probes.addProbeToCollection(PointsProbe(E, Y, std::vector<std::vector<double>>{ {0.3}, { 0.1 } }));

	Sources sources;
	sources.addSourceToVector(Source(model, E, Y, 1.0, 0.5, Vector({ 0.2 })));
	sources.addSourceToVector(Source(model, H, Z, 1.0, 0.5, Vector({ 0.2 })));

	maxwell::Solver solver(
		model,
		probes,
		sources,
		buildDefaultSolverOpts(1.5));

	auto eOld = solver.getFieldInDirection(E, Y);

	double reflectCoeff =
		(model.getAttToMat().at(2).getImpedance() - model.getAttToMat().at(1).getImpedance()) /
		(model.getAttToMat().at(2).getImpedance() + model.getAttToMat().at(1).getImpedance());

	solver.run();

	EXPECT_NEAR(eOld.Max(), 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0, Y), 2e-3);
	EXPECT_NEAR(0.0, 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.45, 0, Y), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0, Y) * reflectCoeff,
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.90, 0, Y), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0, Y) * reflectCoeff,
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.10, 1, Y), 2e-3);
	EXPECT_NEAR(0.0, 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.30, 1, Y), 2e-3);
}
TEST_F(TestMaxwellSolver, twoDimensional_Periodic) //TODO ADD ENERGY CHECK
{
	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
	Vector periodic({ 0.0, 1.0 });
	std::vector<Vector> trans;
	trans.push_back(periodic);
	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D,mesh2D.CreatePeriodicVertexMapping(trans));
	
	Model model = Model(mesh2DPer, AttributeToMaterial(), AttributeToBoundary());
	
	Probes probes;
	probes.addProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

	Sources sources;
	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.2, 0.0 })));

	maxwell::Solver solver(model, probes, sources, buildDefaultSolverOpts(1.0));

	solver.run();

}
TEST_F(TestMaxwellSolver, twoDimensional_centered_NC_MESH) //TODO ADD ENERGY CHECK
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
	probes.addProbeToCollection(ExporterProbe());
	probes.vis_steps = 20;

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
TEST_F(TestMaxwellSolver, twoDimensional_centered_AMR_MESH)
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

//TEST_F(TestMaxwellSolver, DISABLED_twoDimensionalResonantBox)
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