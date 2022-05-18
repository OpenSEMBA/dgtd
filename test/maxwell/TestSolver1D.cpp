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

	AttributeToMaterial buildAttToMatVec(
		const std::vector<Attribute>& attVec, 
		const std::vector<Material>& matVec)
	{
		AttributeToMaterial res;
		for (int i = 0; i < attVec.size(); i++) {
			res.emplace(attVec[i], matVec[i]);
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

class TestMaxwellSolver1D : public ::testing::Test {
protected:

	Model buildOneDimFourMatModel(
		const int meshIntervals) {

		std::vector<Attribute> attArrMultiple = std::vector<Attribute>({ 1, 2, 3, 4 });
		Material mat11 = Material(1.0, 1.0); Material mat12 = Material(1.0, 2.0);
		Material mat21 = Material(2.0, 1.0); Material mat22 = Material(2.0, 2.0);
		std::vector<Material> matArrMultiple;
		matArrMultiple.push_back(mat11); matArrMultiple.push_back(mat12);
		matArrMultiple.push_back(mat21); matArrMultiple.push_back(mat22);
		AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatVec(attArrMultiple, matArrMultiple);
		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), attToMatVec);
	}

	Model buildOneDimOneMatModel(
		const int meshIntervals = 51) {

		std::vector<Attribute> attArrSingle = std::vector<Attribute>({ 1 });
		Material mat11 = Material(1.0, 1.0);
		std::vector<Material> matArrSimple = std::vector<Material>({ mat11 });
		AttributeToMaterial attToMatVec = HelperFunctions::buildAttToMatVec(attArrSingle, matArrSimple);
		return Model(Mesh::MakeCartesian1D(meshIntervals, 1.0), attToMatVec);
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

TEST_F(TestMaxwellSolver1D, checkTwoAttributeMesh)
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
TEST_F(TestMaxwellSolver1D, oneDimensional_centered)
{
	//DEPRECATED INTRO // REWRITE
	
	/*The purpose of this test is to check the run() function for the Solver1D class
	and test the different available options.

	First, dimensional variables are declared and a mesh is constructed, along with the declaration
	of different useful variables.

	Then, a Solver1D object is constructed using said mesh and options, the bounding box for its mesh
	is extracted and an initial condition is applied to one of its variables. (GridFunction Ez_)

	Lastly, the run() function is called.*/

	maxwell::Solver::Options solverOpts;
	
	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.fluxType = FluxType::Centered;
	solverOpts.t_final = 2.0;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 50;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver1D::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver1D::buildOneDimOneMatModel(), probes, 
						sources, solverOpts);
	
	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver1D, oneDimensional_upwind_PEC)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.t_final = 2.0;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 50;

	Sources sources;
	sources.addSourceToVector(TestMaxwellSolver1D::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver1D::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	double error = eOld.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver1D, oneDimensional_upwind_PMC)
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
	sources.addSourceToVector(TestMaxwellSolver1D::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver1D::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	double error = hOld.DistanceTo(hNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver1D, oneDimensional_upwind_SMA)
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
	sources.addSourceToVector(TestMaxwellSolver1D::buildSourceOneDimOneMat());

	maxwell::Solver solver(TestMaxwellSolver1D::buildOneDimOneMatModel(), probes,
		sources, solverOpts);

	GridFunction eOld = solver.getFieldInDirection(E, X);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, X);
	Vector zero = eNew;
	zero = 0.0;
	double error = zero.DistanceTo(eNew);
	EXPECT_NEAR(0.0, error, 2e-3);

}

TEST_F(TestMaxwellSolver1D, TwoSourceWaveTravelsToTheRight_SMA)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.evolutionOperatorOptions.bdrCond = BdrCond::SMA;
	solverOpts.t_final = 0.7;
	solverOpts.dt = 1e-3;

	Probes probes;
	//probes.paraview = true;
	probes.vis_steps = 5;
	probes.extractDataAtPoints = true;
	DenseMatrix pointMat(1, 2);
	pointMat.Elem(0, 0) = 0.5;
	pointMat.Elem(0, 1) = 0.8;
	FieldType fieldToExtract = E;
	Direction directionToExtract = Y;
	Probe probe(fieldToExtract, directionToExtract, pointMat);
	probes.addProbeToVector(probe);

	double spread = 2.0;
	double coeff = 1.0;
	double dev = 0.0;
	Direction d = Y;
	FieldType ft = E;
	Source EYFieldSource = TestMaxwellSolver1D::buildSourceOneDimOneMat();
	Source HZFieldSource = TestMaxwellSolver1D::buildSourceOneDimOneMat(51, spread, coeff, dev, Z, H);
	Sources sources;
	sources.addSourceToVector(EYFieldSource);
	sources.addSourceToVector(HZFieldSource);

	maxwell::Solver solver(TestMaxwellSolver1D::buildOneDimOneMatModel(), probes,
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

TEST_F(TestMaxwellSolver1D, TwoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	maxwell::Solver::Options solverOpts;

	solverOpts.evolutionOperatorOptions = FiniteElementEvolutionNoCond::Options();
	solverOpts.t_final = 1.5;
	solverOpts.dt = 1e-3;

	Probes probes;
	probes.paraview = true;
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
	std::vector<Attribute> attVec = std::vector<Attribute>({ 1, 2 });
	Model model = Model(mesh1D, HelperFunctions::buildAttToMatVec(attVec, matVec));

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