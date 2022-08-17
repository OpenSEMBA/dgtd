#include "gtest/gtest.h"
#include "AnalyticalFunctions1D.h"

#include "maxwell/Solver.h"
#include "GlobalFunctions.h"

using Interval = std::pair<double, double>;

using namespace AnalyticalFunctions1D;
using namespace maxwell;
using namespace mfem;

class TestSolver1D : public ::testing::Test {
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
		res.addPointsProbeToCollection(PointsProbe(fToExtract, dirToExtract,
			std::vector<std::vector<double>>{{0.0},{0.5},{1.0}}));
		return res;
	}

	maxwell::Solver::Options buildDefaultSolverOpts(const double tFinal = 2.0)
	{
		maxwell::Solver::Options res;

		res.evolutionOperatorOptions = FiniteElementEvolution::Options();
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
		const int denseMatPointByOrder)
	{
		auto itpos = findTimeId(probe.getConstFieldMovie(), timeToFind, 1e-6);
		if (itpos == probe.getConstFieldMovie().end()) {
			throw std::exception("Time value has not been found within the specified tolerance.");
		}
		auto FieldValueForTimeAtPoint = itpos->second.at(denseMatPointByOrder).at(probe.getDirection());

		return FieldValueForTimeAtPoint;
	}

	void setAttributeIntervalMesh1D(
		const std::map<Attribute, Interval>& attToInterval,
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
				mesh.SetAttribute((int)i, (int)kv.first);
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
};

TEST_F(TestSolver1D, centered)
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
TEST_F(TestSolver1D, centered_energy)
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
TEST_F(TestSolver1D, upwind_PEC_EX)
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
TEST_F(TestSolver1D, upwind_PEC_EY)
{
	Model model = buildOneDimOneMatModel();

	auto probes = buildProbesWithDefaultPointsProbe(E, Y);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts());

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestSolver1D, upwind_PEC_EZ)
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

	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0,        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(  eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestSolver1D, upwind_PMC_HX)
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
TEST_F(TestSolver1D, upwind_PMC_HY)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::PMC, BdrCond::PMC);

	auto probes = buildProbesWithDefaultPointsProbe(H, Y);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, H, Y),
		buildDefaultSolverOpts()
	);

	GridFunction hOld = solver.getFieldInDirection(H, Y);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Y);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestSolver1D, upwind_PMC_HZ)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::PMC, BdrCond::PMC);

	auto probes = buildProbesWithDefaultPointsProbe(H, Z);
	probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, H, Z),
		buildDefaultSolverOpts());

	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(     0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);
																					  
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(hOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestSolver1D, upwind_SMA_EX)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	maxwell::Solver solver(
		model,
		buildProbesWithDefaultPointsProbe(E, X),
		buildSourcesWithDefaultSource(model, E, X),
		buildDefaultSolverOpts(0.2));

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}
TEST_F(TestSolver1D, DISABLED_upwind_SMA_EY)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	auto probes = buildProbesWithDefaultPointsProbe(E, Y);
	auto probeZ = PointsProbe(E, Z, std::vector<std::vector<double>>({ {0.0},{0.5},{1.0} }));
	auto probeHY = PointsProbe(H, Y, std::vector<std::vector<double>>({ {0.0},{0.5},{1.0} }));
	probes.addPointsProbeToCollection(probeZ);
	probes.addPointsProbeToCollection(probeHY);
	probes.addExporterProbeToCollection(ExporterProbe());
	

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

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.0, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(1), 1.0, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.0, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.0, 1), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(2), 1.0, 2), 2e-3);

}
TEST_F(TestSolver1D, DISABLED_upwind_SMA_EZ)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	auto probes = buildProbesWithDefaultPointsProbe(E, Z);
	//probes.addExporterProbeToCollection(ExporterProbe());

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Z),
		buildDefaultSolverOpts(1.0));

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_GE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_GE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2));

}
TEST_F(TestSolver1D, DISABLED_strong_flux_PEC_EY)
{
	Mesh mesh = Mesh::MakeCartesian1D(51);
	Model model = Model(mesh, AttributeToMaterial(), AttributeToBoundary());

	maxwell::Solver::Options opts;
	opts.evolutionOperatorOptions = FiniteElementEvolution::Options();
	opts.evolutionOperatorOptions.disForm = DisForm::Strong;
	opts.t_final = 0.5;

	Probes probes = buildProbesWithDefaultPointsProbe(E, Y);
	probes.addExporterProbeToCollection(ExporterProbe());
	probes.vis_steps = 5;

	maxwell::Solver solver(
		model,
		probes,
		buildSourcesWithDefaultSource(model, E, Y),
		opts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestSolver1D, DISABLED_weak_strong_flux_comparison)
{
	Model model = buildOneDimOneMatModel();
	
	maxwell::Solver::Options optsWeak;
	optsWeak.evolutionOperatorOptions = FiniteElementEvolution::Options();

	maxwell::Solver solverWeak(
		model,
		buildProbesWithDefaultPointsProbe(E, Y),
		buildSourcesWithDefaultSource(model, E, Y),
		optsWeak);

	maxwell::Solver::Options optsStrong;
	optsStrong.evolutionOperatorOptions = FiniteElementEvolution::Options();
	optsStrong.evolutionOperatorOptions.disForm = DisForm::Strong;

	maxwell::Solver solverStrong(
		model,
		buildProbesWithDefaultPointsProbe(E, Y),
		buildSourcesWithDefaultSource(model, E, Y),
		optsStrong);

	GridFunction eOldWk = solverWeak.getFieldInDirection(E, Y);
	GridFunction eOldSt = solverStrong.getFieldInDirection(E, Y);
	solverWeak.run();
	solverStrong.run();
	GridFunction eNewWk = solverWeak.getFieldInDirection(E, Y);
	GridFunction eNewSt = solverStrong.getFieldInDirection(E, Y);

	double errorOld = eOldWk.DistanceTo(eOldSt);
	double errorNew = eNewWk.DistanceTo(eNewSt);
	EXPECT_NEAR(0.0, errorOld, 2e-3);
	EXPECT_NEAR(0.0, errorNew, 2e-3);

}
TEST_F(TestSolver1D, twoSourceWaveTravelsToTheRight_SMA)
{
	Model model = buildOneDimOneMatModel(51, BdrCond::SMA, BdrCond::SMA);

	Probes probes;
	probes.addPointsProbeToCollection(PointsProbe(E, Y, std::vector<std::vector<double>>{ {0.5}, { 0.8 } }));
	//probes.addExporterProbeToCollection(ExporterProbe());

	Sources sources;
	sources.addSourceToVector(Source(model, E, Y, 2.0, 1.0, Vector({ 0.0 })));
	sources.addSourceToVector(Source(model, H, Z, 2.0, 1.0, Vector({ 0.0 })));

	maxwell::Solver solver(
		model,
		probes,
		sources,
		buildDefaultSolverOpts(0.7));

	solver.run();

	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.3, 1),
				getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0),
				2e-3);

}
TEST_F(TestSolver1D, twoSourceWaveTwoMaterialsReflection_SMA_PEC)
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
	probes.addPointsProbeToCollection(PointsProbe(E, Y, std::vector<std::vector<double>>{ {0.3}, { 0.1 } }));
	//probes.addExporterProbeToCollection(ExporterProbe());

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
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0), 2e-3);
	EXPECT_NEAR(0.0, 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.45, 0), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.90, 0), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.10, 1), 2e-3);
	EXPECT_NEAR(0.0, 
		getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.30, 1), 2e-3);
}
TEST_F(TestSolver1D, fluxOperator_O2)
{

	Model model = buildOneDimOneMatModel(3);

	maxwell::Solver solver(
		model,
		Probes(),
		buildSourcesWithDefaultSource(model, E, Y),
		buildDefaultSolverOpts(0.1));

	auto MSMat = toEigen(*solver.getFEEvol()
		.getInvMassStiffness(FieldType::E, X).get()->SpMat().ToDenseMatrix());
	auto noDirMat =  toEigen(*solver.getFEEvol()
		.getInvMassNoDirFlux(FieldType::E, FieldType::E).get()->SpMat().ToDenseMatrix());
	auto oneDirMat = toEigen(*solver.getFEEvol()
		.getInvMassOneDirFlux(FieldType::E, FieldType::H, X).get()->SpMat().ToDenseMatrix());
}