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
	static const int defaultNumberOfElements{ 51 };

	Model buildModel(
		const int numberOfElements = defaultNumberOfElements, 
		const BdrCond& bdrL = BdrCond::PEC, 
		const BdrCond& bdrR = BdrCond::PEC) {

		return Model(Mesh::MakeCartesian1D(numberOfElements, 1.0), AttributeToMaterial{}, buildAttrToBdrMap1D(bdrL, bdrR));
	}

	BdrCond buildPerfectBoundary(FieldType f) {
		switch (f) {
		case FieldType::E:
			return BdrCond::PEC;
		case FieldType::H:
			return BdrCond::PMC;
		default:
			throw std::runtime_error("Invalid Field type");
		}
	}

	Sources buildGaussianInitialField(
		const FieldType& ft = E,
		const Direction& d = X,
		const double spread = 2.0,
		const double coeff = 1.0,
		const Vector& center = Vector({ 0.5 })) 
	{
		return { GaussianInitialField(ft, d, spread, coeff, center) };
	}

	Sources buildRightTravelingWaveInitialField(const Vector& center)
	{
		return {
			GaussianInitialField(E, Y, 2.0, 1.0, center),
			GaussianInitialField(H, Z, 2.0, 1.0, center)
		};
	}

	PointsProbe buildPointsProbe(
		const FieldType& fToExtract = E, 
		const Direction& dirToExtract = X)
	{
		return { fToExtract, dirToExtract, Points{ {0.0},{0.5},{1.0} } };
	}

	Probes buildProbes(
		const FieldType& f = E,
		const Direction& d = X)
	{
		Probes r{ {buildPointsProbe(f, d)} };
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
		auto itpos = findTimeId(probe.getFieldMovie(), timeToFind, 1e-6);
		if (itpos == probe.getFieldMovie().end()) {
			throw std::exception("Time value has not been found within the specified tolerance.");
		}
		auto FieldValueForTimeAtPoint = itpos->second.at(denseMatPointByOrder);

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

	double getEnergy(const GridFunction& e, const GridFunction& h) 
	{
		return pow(e.Norml2(), 2.0) + pow(h.Norml2(), 2.0);
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
	maxwell::Solver solver{
		buildModel(),
		{ {buildPointsProbe()}, { ExporterProbe{"Centered"} }},
		buildGaussianInitialField(),
		SolverOptions{}.setCentered()
	};
	
	GridFunction eOld{ solver.getFieldInDirection(E, Y) };
	GridFunction hOld = solver.getFieldInDirection(H, Z);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);
	GridFunction hNew = solver.getFieldInDirection(H, Z);

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 2e-3);
	EXPECT_GE(getEnergy(eOld, hOld), getEnergy(eNew, hNew));
}
TEST_F(TestSolver1D, upwind_perfect_boundary_EH_XYZ)
{
	for (const auto& f : { E, H }) {
		for (const auto& x : { X, Y, Z }) {
			
			maxwell::Solver solver{
				buildModel(defaultNumberOfElements, buildPerfectBoundary(f), buildPerfectBoundary(f)),
				{{buildPointsProbe(f, x)} },
				buildGaussianInitialField(f, x),
				SolverOptions()
			};

			GridFunction fOld = solver.getFieldInDirection(f, x);
			solver.run();
			GridFunction fNew = solver.getFieldInDirection(f, x);

			EXPECT_NEAR(0.0, fOld.DistanceTo(fNew), 2e-3);

			const auto pp{ *solver.getPointsProbe(0) };
			EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(pp, 0.5, 0), 2e-3);
			EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(pp, 0.5, 2), 2e-3);
			EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(pp, 1.5, 0), 2e-3);
			EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(pp, 1.5, 2), 2e-3);

			EXPECT_NE(fOld.Max(), getBoundaryFieldValueAtTime(pp, 0.5, 1));
			EXPECT_NE(fOld.Max(), getBoundaryFieldValueAtTime(pp, 1.5, 1));
		}
	}
}
TEST_F(TestSolver1D, upwind_SMA_E_XYZ)
{
	for (const auto& x : { X, Y, Z }) {
		maxwell::Solver solver(
			buildModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
			{ {buildPointsProbe(E, x)} },
			buildGaussianInitialField(E, x),
			SolverOptions{}
		);

		GridFunction eOld = solver.getFieldInDirection(E, x);
		solver.run();
		GridFunction eNew = solver.getFieldInDirection(E, x);

		double error = eOld.DistanceTo(eNew);
		EXPECT_NEAR(0.0, error, 2e-3);
	}
}
TEST_F(TestSolver1D, wave_travelingToTheRight_SMA)
{

	maxwell::Solver solver{
		buildModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
		{ {PointsProbe{E, Y, Points{ {0.5}, { 0.8 } }}} },
		buildRightTravelingWaveInitialField(Vector({ 0.5 })),
		SolverOptions{}.setFinalTime(0.7)
	};

	solver.run();

	EXPECT_NEAR(getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.3, 1),
				getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.0, 0),
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

	Probes probes{
		{ PointsProbe(E, Y, Points{ {0.3}, { 0.1 } }) }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildRightTravelingWaveInitialField(Vector({ 0.5 })),
		SolverOptions{}.setFinalTime(1.5)
	};

	auto eOld = solver.getFieldInDirection(E, Y);

	double reflectCoeff =
		(model.getAttToMat().at(2).getImpedance() - model.getAttToMat().at(1).getImpedance()) /
		(model.getAttToMat().at(2).getImpedance() + model.getAttToMat().at(1).getImpedance());

	solver.run();

	EXPECT_NEAR(eOld.Max(), getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.0, 0), 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.45, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 1.30, 1), 2e-3);
	
	EXPECT_NEAR(getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		        getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.90, 0), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		        getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 1.10, 1), 2e-3);
}

TEST_F(TestSolver1D, DISABLED_strong_flux_PEC_EY)
{
	SolverOptions opts;
	opts.evolutionOperatorOptions.disForm = DisForm::Strong;
	opts.t_final = 0.5;

	maxwell::Solver solver{
		buildModel(),
		{{buildPointsProbe(E, Y)}},
		buildGaussianInitialField(E, Y),
		opts
	};

	GridFunction eOld = solver.getFieldInDirection(E, Y);
	solver.run();
	GridFunction eNew = solver.getFieldInDirection(E, Y);

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.5, 2), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 1.5, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 1.5, 2), 2e-3);

	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 0.5, 1));
	EXPECT_NE(eOld.Max(), getBoundaryFieldValueAtTime(*solver.getPointsProbe(0), 1.5, 1));

}
TEST_F(TestSolver1D, DISABLED_weak_strong_flux_comparison)
{
	auto model{ buildModel() };
	Probes probes{ {buildPointsProbe(E, Y)} };
	auto source{ buildGaussianInitialField(E, Y) };

	maxwell::Solver solverWeak{ 
		model, probes, source,
		SolverOptions{}
	};

	maxwell::Solver solverStrong{
		model, probes, source,
		SolverOptions{}.setStrongForm()
	};

	GridFunction eOldWk = solverWeak.getFieldInDirection(E, Y);
	solverWeak.run();
	GridFunction eNewWk = solverWeak.getFieldInDirection(E, Y);
	
	GridFunction eOldSt = solverStrong.getFieldInDirection(E, Y);
	solverStrong.run();
	GridFunction eNewSt = solverStrong.getFieldInDirection(E, Y);

	EXPECT_NEAR(0.0, eOldWk.DistanceTo(eOldSt), 2e-3);
	EXPECT_NEAR(0.0, eNewWk.DistanceTo(eNewSt), 2e-3);
}
TEST_F(TestSolver1D, DISABLED_fluxOperator_O2)
{
	maxwell::Solver solver{
		buildModel(3),
		Probes(),
		buildGaussianInitialField(E, Y),
		SolverOptions().setFinalTime(0.1)
	};

	auto MSMat = toEigen(*solver.getFEEvol()
		.getInvMassStiffness(FieldType::E, X).get()->SpMat().ToDenseMatrix());
	auto noDirMat =  toEigen(*solver.getFEEvol()
		.getInvMassNoDirFlux(FieldType::E, FieldType::E).get()->SpMat().ToDenseMatrix());
	auto oneDirMat = toEigen(*solver.getFEEvol()
		.getInvMassOneDirFlux(FieldType::E, FieldType::H, X).get()->SpMat().ToDenseMatrix());
}