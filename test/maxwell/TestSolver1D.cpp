#include "gtest/gtest.h"
#include "SourceFixtures.h"
#include "GlobalFunctions.h"

#include "maxwell/Solver.h"

using Interval = std::pair<double, double>;

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

using Solver = maxwell::Solver;

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
		Probes r{ { buildPointsProbe(f, d)} };
		return r;
	}

	Probes buildExportProbes()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
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

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}
};

TEST_F(TestSolver1D, box_pec_centered_flux)
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
		buildExportProbes(),
		buildGaussianInitialField(E, Y),
		SolverOptions{}
			.setTimeStep(2.5e-3)
			.setCentered()
	};
	
	GridFunction eOld{ solver.getFields().E[Y] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Y] };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

TEST_F(TestSolver1D, box_pmc_centered_flux)
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
		buildModel(defaultNumberOfElements, BdrCond::PMC,BdrCond::PMC),
		buildExportProbes(),
		buildGaussianInitialField(H, Z),
		SolverOptions{}
			.setTimeStep(2.5e-3)
			.setCentered()
	};

	GridFunction hOld{ solver.getFields().H[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction hNew{ solver.getFields().H[Z] };

	EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}
TEST_F(TestSolver1D, box_pec_upwind_flux)
{
	maxwell::Solver solver{
		buildModel(),
		buildExportProbes(),
		buildGaussianInitialField(E, Y),
		SolverOptions{}
			.setTimeStep(1e-3)
	};

	GridFunction eOld{ solver.getFields().E[Y] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Y] };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-2);
}

TEST_F(TestSolver1D, box_pmc_upwind_flux)
{
	maxwell::Solver solver{
		buildModel(defaultNumberOfElements, BdrCond::PMC,BdrCond::PMC),
		buildExportProbes(),
		buildGaussianInitialField(H, Z),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setOrder(2)
	};

	GridFunction hOld{ solver.getFields().H[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction hNew{ solver.getFields().H[Z] };

	EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-2);
}
TEST_F(TestSolver1D, box_SMA)
{
	maxwell::Solver solver(
		buildModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
		buildExportProbes(),
		buildGaussianInitialField(E, Y),
		SolverOptions{}
			.setTimeStep(2.5e-3)
			.setFinalTime(1.25)
	);

	solver.run();

	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), 2e-3);
}

//TEST_F(TestSolver1D, DISABLED_upwind_perfect_boundary_EH_XYZ)
//{
//	for (const auto& f : { E, H }) {
//		for (const auto& x : { X, Y, Z }) {
//			const auto y{ (x + 1) % 3 };
//			maxwell::Solver solver{
//				buildModel(defaultNumberOfElements, buildPerfectBoundary(f), buildPerfectBoundary(f)),
//				{{buildPointsProbe(f, y)} },
//				buildGaussianInitialField(f, y),
//				SolverOptions()
//			};
//
//			GridFunction fOld = solver.getFieldInDirection(f, y);
//			solver.run();
//			GridFunction fNew = solver.getFieldInDirection(f, y);
//
//			EXPECT_NEAR(0.0, fOld.DistanceTo(fNew), 2e-3);
//		}
//	}
//}
TEST_F(TestSolver1D, DISABLED_box_upwind_SMA_E_XYZ)
{
	for (const auto& x : { X, Y, Z }) {
		maxwell::Solver solver(
			buildModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
			{ {buildPointsProbe(E, x)} },
			buildGaussianInitialField(E, x),
			SolverOptions{}
		);

		GridFunction eOld{ solver.getFields().E[x] };
		solver.run();
		GridFunction eNew{ solver.getFields().E[x] };

		double error = eOld.DistanceTo(eNew);
		EXPECT_NEAR(0.0, error, 2e-3);
	}
}

TEST_F(TestSolver1D, DISABLED_twoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	Mesh mesh1D = Mesh::MakeCartesian1D(101);
	setAttributeIntervalMesh1D({ { 2, std::make_pair(0.76, 1.0) } }, mesh1D);
	
	Material mat1{1.0, 1.0};
	Material mat2{2.0, 1.0};

	Model model = Model(
		mesh1D, 
		{ {1, mat1}, {2, mat2} },
		{ {1, BdrCond::SMA}, {2, BdrCond::PEC} }
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

	auto eOld = solver.getFields().E[Y];

	double reflectCoeff =
		(mat2.getImpedance() - mat1.getImpedance()) /
		(mat2.getImpedance() + mat1.getImpedance());

	solver.run();

	EXPECT_NEAR(eOld.Max(), getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0), 2e-3);

	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.45, 0), 2e-3);
	EXPECT_NEAR(0.0, getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.30, 1), 2e-3);
	
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.90, 0), 2e-3);
	EXPECT_NEAR(getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 0.0, 0) * reflectCoeff,
		        getBoundaryFieldValueAtTime(solver.getPointsProbe(0), 1.10, 1), 2e-3);
}