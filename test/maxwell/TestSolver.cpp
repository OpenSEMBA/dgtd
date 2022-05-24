#include "gtest/gtest.h"
#include <math.h>

#include "maxwell/Solver.h"

using namespace maxwell;

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

	Mesh makeTwoAttributeCartesianMesh1D(
		const int& refTimes = 0)
	{
		Mesh res = Mesh::MakeCartesian1D(2);
		res.SetAttribute(0, 1);
		res.SetAttribute(1, 2);

		for (int i = 0; i < refTimes; i++) {
			res.UniformRefinement();
		}

		return res;
	}

	void setAttributeIntervalMesh(
		const int& attVal, 
		const Array<int>& elIndexes, 
		Mesh& mesh)
	{
		if (elIndexes.begin() > elIndexes.end()) {
			throw std::exception("Lower Index bigger than Higher Index.");
		}
		if (elIndexes[1] > mesh.GetNE()) {
			throw std::exception("Declared element index bigger than Mesh Number of Elements.");
		}
		for (int i = elIndexes[0]; i <= elIndexes[1]; i++) {
			mesh.SetAttribute(i, attVal);
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

	AttributeToMaterial buildAttToMatMap(
		const std::vector<Attribute>& attVec, 
		const std::vector<Material>& matVec)
	{

		if (attVec.size() != matVec.size()) {
			throw std::exception("attVec and matVec must have the same size.");
		}
		AttributeToMaterial res;
		for (int i = 0; i < attVec.size(); i++) {
			res.emplace(attVec[i], matVec[i]);
		}
		return res;
	}

	AttributeToBoundary buildAttToBdrMap(
		const std::vector<Attribute>& attVec,
		const std::vector<BdrCond>& bdrVec)
	{
		AttributeToBoundary res;
		for (int i = 0; i < attVec.size(); i++) {
			res.emplace(attVec[i], bdrVec[i]);
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

class TestMaxwellSolver : public ::testing::Test {
protected:

	Model buildOneDimFourMatModel(
		const int meshIntervals) {

		std::vector<Attribute> attArrMultiple = std::vector<Attribute>({ 1, 2, 3, 4 });
		Material mat11 = Material(1.0, 1.0); Material mat12 = Material(1.0, 2.0);
		Material mat21 = Material(2.0, 1.0); Material mat22 = Material(2.0, 2.0);
		std::vector<Material> matArrMultiple;
		matArrMultiple.push_back(mat11); matArrMultiple.push_back(mat12);
		matArrMultiple.push_back(mat21); matArrMultiple.push_back(mat22);
		AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatMap(attArrMultiple, matArrMultiple);
		AttributeToBoundary attToBdrVec;
		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), attToMatVec, attToBdrVec);
	}

	Model buildOneDimOneMatModel(
		const int meshIntervals = 51) {

		std::vector<Attribute> attArrSingle = std::vector<Attribute>({ 1 });
		Material mat11 = Material(1.0, 1.0);
		std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
		AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatMap(attArrSingle, matArrSimple);
		AttributeToBoundary attToBdrVec;
		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), attToMatVec, attToBdrVec);
	}

	Source buildSourceOneDimOneMat(
		const int meshIntervals = 51, 
		const double spread = 2.0, 
		const double coeff = 1.0, 
		const double dev = 0.0,
		const Direction& d = Y, 
		const FieldType& ft = E){

		return Source(buildOneDimOneMatModel(meshIntervals), spread, coeff, dev, d, ft);
	}

};

TEST_F(TestMaxwellSolver, checkTwoAttributeMesh)
{
	/*The purpose of this test is to check the makeTwoAttributeCartesianMesh1D(const int& refTimes) 
	function.

	First, an integer is declared for the number of times we wish to refine the mesh, then a mesh is 
	constructed with two elements, left and right hand sides, setting the following attributes.

	|------LHS------|------RHS------|

	|##ATTRIBUTE 1##|##ATTRIBUTE 2##|

	Once the mesh is refined, it is returned, then we compare if the expected number of elements is
	true for the actual elements in the mesh.

	Then, we consider how the mesh will perform its uniform refinement, and we declare that the 
	LHS elements with Attribute one will be Even index elements (starting at 0), and the RHS
	elements with Attribute 2 will be Uneven index elements (starting at 1).*/
	
	const int refTimes = 3;
	Mesh mesh = HelperFunctions::makeTwoAttributeCartesianMesh1D(refTimes);

	EXPECT_EQ(pow(2,refTimes + 1), mesh.GetNE());
	for (int i = 0; i < mesh.GetNE(); i++) {
		if (i % 2 == 0) {
			EXPECT_EQ(1, mesh.GetAttribute(i));
		}
		else {
			EXPECT_EQ(2, mesh.GetAttribute(i));
		}
	}
}

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

	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.t_final = 2.0;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 50;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver::buildOneDimOneMatModel(), probes, 
						sources, solverOpts);
	
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

	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.t_final = 1.999;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 50;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	EXPECT_GE(eOld.Norml2() + hOld.Norml2(), eNew.Norml2() + hNew.Norml2());

}
TEST_F(TestMaxwellSolver, oneDimensional_upwind_PEC)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.t_final = 2.0;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 50;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver, oneDimensional_upwind_PMC)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::PMC;
	solverOpts.t_final = 1.0;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 5;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver, oneDimensional_upwind_SMA)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::SMA;
	solverOpts.t_final = 1.0;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 5;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver, twoSourceWaveTravelsToTheRight_SMA)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::SMA;
	solverOpts.t_final = 0.7;
	solverOpts.dt = 1e-3;

	Probes probes;
	probes.paraview = true;
	probes.vis_steps = 5;
	//probes.extractDataAtPoints = true;
	DenseMatrix pointMat(1, 2);
	pointMat.Elem(0, 0) = 0.5;
	pointMat.Elem(0, 1) = 0.8;
	FieldType fieldToExtract = E;
	Direction directionToExtract = Y;
	Probe probe(fieldToExtract, directionToExtract, pointMat);
	probes.addProbeToVector(probe);


	std::vector<Material> matVec;
	matVec.push_back(Material(1.0, 1.0));
	std::vector<Attribute> attVec = std::vector<Attribute>({1});
	std::vector<BdrCond> bdrVec;
	bdrVec.push_back(BdrCond::SMA);
	bdrVec.push_back(BdrCond::SMA);
	std::vector<Attribute> bdrAttVec = std::vector<Attribute>({ 1, 2 });
	Model model = Model(Mesh::MakeCartesian1D(51), HelperFunctions::buildAttToMatMap(attVec, matVec), HelperFunctions::buildAttToBdrMap(bdrAttVec, bdrVec));

	double spread = 2.0;
	double coeff = 1.0;
	double dev = 0.0;
	Direction d = Y;
	FieldType ft = E;
	Source EYFieldSource = Source(model, spread, coeff, dev, d, ft);
	Source HZFieldSource = Source(model, spread, coeff, dev, Z, H);
	Sources sources;
	sources.addSourceToVector(EYFieldSource);
	sources.addSourceToVector(HZFieldSource);

	maxwell::Solver solver(model, probes,
		sources, solverOpts);

	///////////////////

	solver.run();
	Probe probeEY = solver.getProbe(0);

	///////////////////

	auto it = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 0.0, 1e-6);
	if (it == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValForFirstPos = it->second.at(0).at(Y);

	auto it2 = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 0.3, 1e-6);
	if (it2 == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValForSecondPos = it2->second.at(1).at(Y);
	
	EXPECT_NEAR(EYValForSecondPos, EYValForFirstPos, 2e-3);

}

TEST_F(TestMaxwellSolver, twoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.t_final = 1.5;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 10;
	probes.extractDataAtPoints = true;
	DenseMatrix pointMat(1, 2);
	pointMat.Elem(0, 0) = 0.3;
	pointMat.Elem(0, 1) = 0.1;
	FieldType fieldToExtract = E;
	Direction directionToExtract = Y;
	Probe probe(fieldToExtract, directionToExtract, pointMat);
	probes.addProbeToVector(probe);

	int meshIntervals = 101;
	Mesh mesh1D = Mesh::MakeCartesian1D(meshIntervals);
	DenseMatrix changeAttMat(1,2);
	changeAttMat.Elem(0, 0) = 0.76;
	changeAttMat.Elem(0, 1) = 1.0;
	Array<int> elemID;
	Array<IntegrationPoint> integPoint;
	mesh1D.FindPoints(changeAttMat, elemID, integPoint);

	HelperFunctions::setAttributeIntervalMesh(2,elemID,mesh1D);

	std::vector<Material> matVec;
	matVec.push_back(Material(1.0, 1.0));
	matVec.push_back(Material(2.0, 1.0));
	std::vector<Attribute> matAttVec = std::vector<Attribute>({ 1, 2 });
	std::vector<BdrCond> bdrCondVec;
	bdrCondVec.push_back(BdrCond::SMA);
	bdrCondVec.push_back(BdrCond::PEC);	
	std::vector<Attribute> bdrAttVec = std::vector<Attribute>({ 1, 2 });
	Model model = Model(mesh1D, HelperFunctions::buildAttToMatMap(matAttVec, matVec), HelperFunctions::buildAttToBdrMap(bdrAttVec,bdrCondVec));

	double spread = 1.0;
	double coeff = 0.5;
	double dev = 0.2;
	Direction d = Y;
	FieldType ft = E;

	Source EYFieldSource = Source(model, spread, coeff, dev, d, ft);
	Source HZFieldSource = Source(model, spread, coeff, dev, Z, H);
	Sources sources;
	sources.addSourceToVector(EYFieldSource);
	sources.addSourceToVector(HZFieldSource);

	maxwell::Solver solver(model, probes,
		sources, solverOpts);

	///////////////////

	solver.run();
	Probe probeEY = solver.getProbe(0);

	///////////////////

	auto it0 = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 0.0, 1e-6);
	if (it0 == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValAtTime0 = it0->second.at(0).at(Y);
	
	auto it45 = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 0.45, 1e-6);
	if (it45 == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValAtTime45 = it45->second.at(0).at(Y);

	auto it90 = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 0.90, 1e-6);
	if (it90 == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValAtTime90 = it90->second.at(0).at(Y);

	auto it110 = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 1.10, 1e-6);
	if (it110 == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValAtTime110 = it110->second.at(1).at(Y);

	auto it130 = HelperFunctions::findTimeId(probeEY.getFieldMovie(), 1.30, 1e-6);
	if (it130 == probeEY.getFieldMovie().end()) {
		GTEST_FATAL_FAILURE_("Time value has not been found within the specified tolerance.");
	}
	auto EYValAtTime130 = it130->second.at(1).at(Y);

	double reflectCoeff =
		(matVec.at(1).getImpedance() - matVec.at(0).getImpedance()) / 
		((matVec.at(1).getImpedance() + matVec.at(0).getImpedance()));

	EXPECT_NEAR(0.0, EYValAtTime45, 2e-3);
	EXPECT_NEAR(EYValAtTime0 * reflectCoeff, EYValAtTime90, 2e-3);
	EXPECT_NEAR(EYValAtTime0 * reflectCoeff, EYValAtTime110, 2e-3);
	EXPECT_NEAR(0.0, EYValAtTime130, 2e-3);

}

TEST_F(TestMaxwellSolver, twoDimensionalResonantBox)
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
	double dev = 0.0;
	Source EXFieldSource = Source(model, spread, coeff, dev, X, E); 
	Sources sources;
	sources.addSourceToVector(EXFieldSource);

	Probes probes;
	probes.paraview = true;
	probes.vis_steps = 100;

	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.t_final = 2.0;
	solverOpts.dt = 1e-4;
	solverOpts.order = 1;

	maxwell::Solver solver(model, probes,
		sources, solverOpts);

	solver.run();

}


TEST_F(TestMaxwellSolver, twoDimensional_centered_NCMESH)
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

	maxwell::Solver::Options solverOpts;
	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.t_final = 2.92;
	solverOpts.dt = 1e-3;

	Probes probes;
	probes.paraview = true;
	probes.vis_steps = 20;

	const char* mesh_file = "star-mixed.mesh";
	Mesh mesh(mesh_file);
	std::vector<Attribute> attArrSingle = std::vector<Attribute>({ 1 });
	Material mat11 = Material(1.0, 1.0);
	std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
	AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatMap(attArrSingle, matArrSimple);
	AttributeToBoundary attToBdrVec;
	Model model = Model(mesh, attToMatVec, attToBdrVec);

	double spread = 2.0;
	double coeff = 20.0;
	double dev = 0.0;
	Source EXFieldSource = Source(model, spread, coeff, dev, Z, E);
	Sources sources;
	sources.addSourceToVector(EXFieldSource);

	maxwell::Solver solver(model, probes,
		sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Z);

	EXPECT_GT(eOld.Max(), eNew.Max());
}
//
//TEST_F(TestMaxwellSolver1D, oneDimensional_two_materials)
//{
//	int nx = 100;
//	mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(nx);
//	HelperFunctions1D::setAttributeIntervalMesh(2, Vector({ 51,100 }), mesh);
//
//	maxwell::Solver1D::Options solverOpts;
//	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
//	solverOpts.t_final = 2.999;
//	
//	Vector eps = Vector({ 1.0, 2.0 });
//	Vector mu = Vector({ 1.0, 1.0 });
//	std::list<Material> matArray;
//	for (int i = 0; i < eps.Size(); i++) { //Check if mesh.attributes.Max() broken
//		Material matAux(eps[i], mu[i]);
//		matArray.push_back(matAux);
//	}
//
//	maxwell::Solver1D solver(solverOpts, mesh);
//	
//	solver.getMesh().GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax);
//	solver.setInitialField(FieldType::E, gaussianFunctionHalfWidth);
//
//	Vector eOld = solver.getField(FieldType::E);
//	solver.run();
//	Vector eNew = solver.getField(FieldType::E);
//
//	double error = eOld.DistanceTo(eNew);
//	EXPECT_NEAR(0.0, error, 2e-3);
//}