#include <gtest/gtest.h>

#include "TestUtils.h"

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"

using Interval = std::pair<double, double>;

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

using Solver = maxwell::Solver;

class Solver1DTest : public ::testing::Test {
protected:
	static const int defaultNumberOfElements{ 51 };

	Model buildStandardModel(
		const int numberOfElements = defaultNumberOfElements, 
		const BdrCond& bdrL = BdrCond::PEC, 
		const BdrCond& bdrR = BdrCond::PEC) {
		auto msh{ Mesh::MakeCartesian1D(numberOfElements, 1.0) };
		return Model{
			msh,
			GeomTagToMaterialInfo(),
			GeomTagToBoundaryInfo(GeomTagToBoundary{{1, bdrL},{2, bdrR}}, GeomTagToInteriorBoundary{})
		};
	}

	void setAttributeOnInterval(
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
				throw std::runtime_error("Lower Index bigger than Higher Index.");
			}
			if (elemID[1] > mesh.GetNE()) {
				throw std::runtime_error("Declared element index bigger than Mesh Number of Elements.");
			}
			for (int i = elemID[0]; i <= elemID[1]; i++) {
				mesh.SetAttribute((int)i, (int)kv.first);
			}
		}
	}

	void expectFieldsAreNearAfterEvolution(maxwell::Solver& solver)
	{
		GridFunction eOld{ solver.getField(E,Y) };
		GridFunction hOld{ solver.getField(H,Z) };

		solver.run();

		GridFunction eNew{ solver.getField(E,Y) };
		GridFunction hNew{ solver.getField(H,Z) };

		EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
		EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	}

};

TEST_F(Solver1DTest, pec_centered)
{
	// This test checks propagation of a wave inside a PEC box. 
	// Final time is set so that a full cycle is completed.

	// auto probes{ buildProbesWithAnExportProbe(50) }; // For DEBUGGING.
	auto probes{ buildProbesEmpty() };
	probes.fieldProbes = {
		FieldProbe{E, Y, {0.0}},
		FieldProbe{H, Z, {0.0}}
	};

	double tolerance{ 1e-3 };

	maxwell::Solver solver{
		buildStandardModel(),
		probes,
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}
			.setCentered()
	};
	
	GridFunction eOld{ solver.getField(E,Y) };
	GridFunction hOld{ solver.getField(H,Z) };

	auto eNormOld{ solver.getFields().getNorml2() };
	
	// Checks fields have been initialized.
	EXPECT_NE(0.0, eNormOld); 
	
	solver.run();
	
	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getField(E,Y) };
	GridFunction hNew{ solver.getField(H,Z) };

	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);

	EXPECT_NEAR(eNormOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary the electric field should be always close to zero...
	for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
		EXPECT_NEAR(0.0, f, tolerance);
	}
	
	// ... and the magnetic field reaches a maximum close to 1.0 
	// (the wave splits in two and doubles at the boundary).
	auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	EXPECT_NEAR(1.5, hMaxFrame.first, 0.01);
	EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);
}

TEST_F(Solver1DTest, pmc_centered)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single Gaussian is declared along Ey.

	Then, the Solver object is constructed using said parts, with its m being one-dimensional.
	The field along Ey is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that after two seconds with PEC boundary conditions, the wave evolves
	back to its initial state within the specified error.*/
	maxwell::Solver solver{
		buildStandardModel(defaultNumberOfElements, BdrCond::PMC, BdrCond::PMC),
		buildProbesEmpty(),
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}.setCentered()
	};

	expectFieldsAreNearAfterEvolution(solver);
}

TEST_F(Solver1DTest, pec_upwind)
{
	maxwell::Solver solver{
		buildStandardModel(),
		buildProbesEmpty(),
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}
	};

	expectFieldsAreNearAfterEvolution(solver);
}

TEST_F(Solver1DTest, pmc_upwind)
{
	maxwell::Solver solver{
		buildStandardModel(defaultNumberOfElements, BdrCond::PMC, BdrCond::PMC),
		buildProbesEmpty(),
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}
	};

	expectFieldsAreNearAfterEvolution(solver);
}

TEST_F(Solver1DTest, sma)
{
	maxwell::Solver solver(
		buildStandardModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
		buildProbesEmpty(),
		buildGaussianInitialField(E, 0.1, Vector({ 0.5 }), unitVec(Y)),
		SolverOptions{}
	);

	EXPECT_NE(0.0, solver.getFields().get(E,Y).Norml2());
	EXPECT_NEAR(0.0, solver.getFields().get(H,Z).Norml2(), 2e-3);

	solver.run();

	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), 2e-3);
}

TEST_F(Solver1DTest, DISABLED_periodic)
{
	auto m{ 
		Mesh::LoadFromFile((mfemMeshes1DFolder() + "periodic-segment.mesh"), 1, 0)
	};

	for (auto i{0}; i < 4; ++i) {
		m.UniformRefinement();
	}

	maxwell::Solver solver{
		Model{ m },
		buildProbesEmpty(),
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}.setFinalTime(1.0)
	};

	GridFunction eOld{ solver.getField(E,Y) };
	GridFunction hOld{ solver.getField(H,Z) };

	solver.run();
	{	
		GridFunction eNew{ solver.getField(E,Y) };
		GridFunction hNew{ solver.getField(H,Z) };
		EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
		EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	}

	solver.setFinalTime(2.0);
	solver.run();
	{	
		GridFunction eNew{ solver.getField(E,Y) };
		GridFunction hNew{ solver.getField(H,Z) };
		EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
		EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	}

}

TEST_F(Solver1DTest, DISABLED_periodic_inhomo)
{
	Mesh m{ Mesh::LoadFromFile(
		(mfemMeshes1DFolder() + "periodic-inhomo-segment.mesh"),1,0) 
	};
	for (auto i{0}; i < 2; ++i) {
		m.UniformRefinement();
	}

	maxwell::Solver solver{
		{m},
		buildProbesEmpty(),
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}
	};

	expectFieldsAreNearAfterEvolution(solver);
}

TEST_F(Solver1DTest, twoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	// Sends a wave through a material interface. 
	// Checks reflection and transmission.
	// Ref: https://en.wikipedia.org/wiki/Reflection_coefficient
	
	auto msh{ Mesh::MakeCartesian1D(100) };

	setAttributeOnInterval({ { 2, std::make_pair(0.50, 1.0) } }, msh);
	
	Material mat1{1.0, 1.0, 0.0};
	Material mat2{4.0, 1.0, 0.0};

	auto probes{ buildProbesWithAnExportProbe(100) };
	probes.fieldProbes = {
		FieldProbe{ E, Y, {0.00} },
		FieldProbe{ E, Y, {0.75} }
	};
	
	maxwell::Solver solver{
		Model{
			msh,
			GeomTagToMaterialInfo(GeomTagToMaterial{ {1, mat1}, {2, mat2} }, GeomTagToBoundaryMaterial{}),
			GeomTagToBoundaryInfo(GeomTagToBoundary{ {1, BdrCond::SMA}, {2, BdrCond::PEC} }, GeomTagToInteriorBoundary{})
		},
		probes,
		buildPlanewaveInitialField(
			Gaussian{ 0.05,  Position(0.0) },
			Source::Position({ 0.25 }),
			Source::Polarization(unitVec(Y)),
			Source::Propagation(unitVec(X))
		),
		SolverOptions{}.setFinalTime(1.0)
	};
		
	solver.run();

	auto expectedReflectCoeff{
		(mat2.getImpedance() - mat1.getImpedance()) /
		(mat2.getImpedance() + mat1.getImpedance())
	};
	auto expectedTransmissionCoeff{ 1 + expectedReflectCoeff };

	auto timeTolerance{ 0.03 };
	auto fieldTolerance{ 0.01 };

	// Checks reflected wave.
	{
		auto frame{ solver.getPointProbe(0).findFrameWithMin() };
		EXPECT_NEAR(0.75, frame.first, timeTolerance);
		EXPECT_NEAR(expectedReflectCoeff, frame.second, fieldTolerance);
	}

	// Checks transmitted wave.
	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		auto expectedTimeOfArrival{ 0.25 + 0.25 / mat2.getSpeedOfWave() };
		EXPECT_NEAR(expectedTimeOfArrival, frame.first, timeTolerance);
		EXPECT_NEAR(expectedTransmissionCoeff, frame.second, fieldTolerance);
	}

}