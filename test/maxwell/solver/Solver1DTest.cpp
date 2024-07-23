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
			GeomTagToMaterial{},
			GeomTagToBoundary{
				{1, bdrL},
				{2, bdrR}
			}
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
	probes.pointProbes = {
		PointProbe{E, Y, {0.0}},
		PointProbe{H, Z, {0.0}}
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
			.setOrder(3)
	);

	EXPECT_NEAR(6.0131571207153787, solver.getFields().getNorml2(), 2e-3);

	solver.run();

	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), 2e-3);
}

TEST_F(Solver1DTest, periodic)
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

TEST_F(Solver1DTest, periodic_inhomo)
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
	
	Material mat1{1.0, 1.0};
	Material mat2{4.0, 1.0};

	auto probes{ buildProbesWithAnExportProbe(100) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.00} },
		PointProbe{ E, Y, {0.75} }
	};
	
	maxwell::Solver solver{
		Model{
			msh,
			{ {1, mat1}, {2, mat2} },
			{ {1, BdrCond::SMA}, {2, BdrCond::PEC} }
		},
		probes,
		buildPlanewaveInitialField(
			Gaussian{ 0.05 },
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
		auto expectedTimeOfArrival{ 0.25 + 0.25 / mat2.getSpeedOfLight() };
		EXPECT_NEAR(expectedTimeOfArrival, frame.first, timeTolerance);
		EXPECT_NEAR(expectedTransmissionCoeff, frame.second, fieldTolerance);
	}

}

TEST_F(Solver1DTest, totalfieldin_intbdr_centered)
{
	auto mesh{ 
		Mesh::LoadFromFile(
			(mfemMeshes1DFolder() + "longlineIntBdr.mesh"), 1, 0
		)
	};
	GeomTagToBoundary attToBdr{ {2, BdrCond::PEC} };
	GeomTagToInteriorConditions attToIntBdr{ {301,BdrCond::TotalFieldIn} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe(50) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.1001} },
		PointProbe{ E, Y, {1.0} },
		PointProbe{ H, Z, {1.0} }
	};
	
	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Y), unitVec(X)),
		SolverOptions{}
			.setCentered()
			.setFinalTime(4.0)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).findFrameWithMax() };
		EXPECT_NEAR(1.5, frame.first, 1e-1);
		EXPECT_NEAR(1.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}
	
	{
		auto frame{ solver.getPointProbe(2).findFrameWithMax() };
		EXPECT_NEAR(2.5, frame.first, 2e-1);
		EXPECT_NEAR(2.0, frame.second, 1e-3);
	}
}

TEST_F(Solver1DTest, totalfieldin_intbdr_submesher_centered)
{
	auto mesh{
		Mesh::LoadFromFile(
			(mfemMeshes1DFolder() + "longlineIntBdr.mesh"), 1, 0
		)
	};
	GeomTagToBoundary attToBdr{ {2, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(20) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.1001} },
		PointProbe{ E, Y, {1.0} },
		PointProbe{ H, Z, {1.0} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Y), unitVec(X)),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(4.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).findFrameWithMax() };
		EXPECT_NEAR(1.5, frame.first, 1e-1);
		EXPECT_NEAR(1.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).findFrameWithMax() };
		EXPECT_NEAR(2.5, frame.first, 2e-1);
		EXPECT_NEAR(2.0, frame.second, 1e-3);
	}
}

TEST_F(Solver1DTest, totalfieldin_intbdr_submesher_upwind)
{
	auto mesh{
		Mesh::LoadFromFile(
			(mfemMeshes1DFolder() + "longlineIntBdr.mesh"), 1, 0
		)
	};
	GeomTagToBoundary attToBdr{ {2, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(20) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.1001} },
		PointProbe{ E, Y, {1.0} },
		PointProbe{ H, Z, {1.0} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Y), unitVec(X)),
		SolverOptions{}
			.setCFL(0.4)
			.setFinalTime(4.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).findFrameWithMax() };
		EXPECT_NEAR(1.5, frame.first, 1e-1);
		EXPECT_NEAR(1.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).findFrameWithMax() };
		EXPECT_NEAR(2.5, frame.first, 2e-1);
		EXPECT_NEAR(2.0, frame.second, 1e-3);
	}
}

TEST_F(Solver1DTest, totalfield_doublebdr_submesher_centered)
{
	auto mesh{
		Mesh::LoadFromFileNoBdrFix(
			(gmshMeshesFolder() + "1D_TFp_line.msh"), 1, 0
		)
	};
	GeomTagToBoundary attToBdr{ {2, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(5) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.5001} },
		PointProbe{ E, Y, {2.0} },
		PointProbe{ H, Z, {2.0} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Y), unitVec(X)),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(4.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).findFrameWithMax() };
		EXPECT_NEAR(2.0, frame.first, 1e-1);
		EXPECT_NEAR(1.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).findFrameWithMax() };
		EXPECT_NEAR(2.5, frame.first, 2e-1);
		EXPECT_NEAR(2.0, frame.second, 1e-3);
	}
}

TEST_F(Solver1DTest, totalfieldinout_intbdr_submesher_centered)
{
	Mesh mesh{	Mesh::LoadFromFileNoBdrFix((mfemMeshes1DFolder() + "LineTFSFInOut.mesh"), 1, 0)};
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(20) };
	probes.pointProbes = {
	PointProbe{ E, Y, {0.1001} },
	PointProbe{ E, Y, {1.0} },
	PointProbe{ H, Z, {0.9} },
	PointProbe{ H, Z, {1.0} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Y), unitVec(X)),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).getFieldMovie() };
		auto expected_t = 1.6;
		for (const auto& [t, f] : frame) {
			if (abs(t - expected_t) <= 1e-2) {
				EXPECT_NEAR(f, 1.0, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).getFieldMovie() };
		auto expected_t = 2.4;
		for (const auto& [t, f] : frame) {
			if (abs(t - expected_t) <= 1e-2) {
				EXPECT_NEAR(f, 1.0, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(3).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}
}

TEST_F(Solver1DTest, totalfieldinout_rtl_intbdr_submesher_centered)
{
	Mesh mesh{ Mesh::LoadFromFileNoBdrFix((mfemMeshes1DFolder() + + "LineTFSFInOut_RtL.mesh"), 1, 0) };
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };
	auto probes{ buildProbesWithAnExportProbe(20) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.9} },
		PointProbe{ E, Y, {0.95} },
		PointProbe{ H, Z, {0.2} },
		PointProbe{ H, Z, {0.95} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 2.4, unitVec(Y), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).getFieldMovie()};
		auto expected_t = 1.5;
		for (const auto& [t, f] : frame) {
			if (abs(t - expected_t) <= 1e-2) {
				EXPECT_NEAR(f, 1.0, 1e-3);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).getFieldMovie()};
		auto expected_t = 2.2;
		for (const auto& [t, f] : frame) {
			if (abs(t - expected_t) <= 1e-2) {
				EXPECT_NEAR(f, -1.0, 1e-3);
			}
		}
	}


	{
		auto frame{ solver.getPointProbe(3).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

}

TEST_F(Solver1DTest, totalfieldinout_rtl_intbdr_submesher_upwind)
{
	Mesh mesh{ Mesh::LoadFromFileNoBdrFix((mfemMeshes1DFolder() + +"LineTFSFInOut_RtL.mesh"), 1, 0) };
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };
	auto probes{ buildProbesWithAnExportProbe(20) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.9} },
		PointProbe{ E, Y, {0.95} },
		PointProbe{ H, Z, {0.2} },
		PointProbe{ H, Z, {0.95} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 2.4, unitVec(Y), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setCFL(0.5)
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).getFieldMovie() };
		auto expected_t = 1.5;
		for (const auto& [t, f] : frame) {
			if (abs(t - expected_t) <= 1e-2) {
				EXPECT_NEAR(f, 1.0, 1e-3);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).getFieldMovie() };
		auto expected_t = 2.2;
		for (const auto& [t, f] : frame) {
			if (abs(t - expected_t) <= 1e-2) {
				EXPECT_NEAR(f, -1.0, 1e-3);
			}
		}
	}


	{
		auto frame{ solver.getPointProbe(3).findFrameWithMax() };
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

}

TEST_F(Solver1DTest, totalfieldin_shortline_intbdr_submesher_centered)
{
	Mesh mesh{
		Mesh::LoadFromFile(
			(mfemMeshes1DFolder() + "lineIntBdr.mesh"), 1, 0
		)
	};
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, attToBdr, GeomTagToInteriorConditions{} };

	auto probes{ buildProbesWithAnExportProbe(1) };
	probes.pointProbes = {
		PointProbe{ E, Y, {0.1001} },
		PointProbe{ E, Y, {2.0} },
		PointProbe{ H, Z, {2.0} },
		PointProbe{ H, Z, {0.1001} }
	};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Y), unitVec(X)),
		SolverOptions{}
			.setCFL(0.1)
			.setCentered()
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	{
		auto frame{ solver.getPointProbe(0).findFrameWithMax() };
		EXPECT_NEAR(3.5, frame.first, 1e-1);
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		EXPECT_NEAR(3.5, frame.first, 1e-1);
		EXPECT_NEAR(1.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(2).findFrameWithMax() };
		EXPECT_NEAR(3.5, frame.first, 1e-1);
		EXPECT_NEAR(1.0, frame.second, 1e-3);
	}

	{
		auto frame{ solver.getPointProbe(3).findFrameWithMax() };
		EXPECT_NEAR(3.5, frame.first, 1e-3);
		EXPECT_NEAR(0.0, frame.second, 1e-3);
	}

}

TEST_F(Solver1DTest, DISABLED_resonant_mode_upwind)
{
	// Resonant mode inside a PEC box. 
	auto probes{ buildProbesWithAnExportProbe() };

	double finalTime{ 1.2 };

	maxwell::Solver solver{
		buildStandardModel(),
		probes,
		buildResonantModeInitialField(E, unitVec(Y), {1}),
		SolverOptions{}
			.setFinalTime(finalTime)
			.setCFL(0.5)
	};

	Vector eOld(solver.getField(E,Y).Size());
	eOld = solver.getField(E,Y);
	
	solver.run();

	GridFunction eNew{ solver.getField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-8);

	//EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	//for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
	//	EXPECT_NEAR(0.0, f, tolerance);
	//}

	//auto hMaxFrame{ solver.getPointProbe(1).findFrameWithMax() };
	//EXPECT_NEAR(1.5, hMaxFrame.first, 0.01);
	//EXPECT_NEAR(1.0, hMaxFrame.second, tolerance);

	EXPECT_TRUE(false);
}

TEST_F(Solver1DTest, pec_centered_spectral)
{
	// This test checks propagation of a wave inside a PEC box. 
	// Final time is set so that a full cycle is completed.
	auto probes{ buildProbesWithAnExportProbe() };
	probes.pointProbes = {
		PointProbe{E, Y, {0.0}},
		PointProbe{H, Z, {0.0}}
	};

	double tolerance{ 1e-3 };

	maxwell::Solver solver{
		buildStandardModel(),
		probes,
		buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
		SolverOptions{}
			.setCentered()
			.setSpectralEO(true)
	};

	GridFunction eOld{ solver.getField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	solver.run();

	// Checks that field is almost the same as initially because the completion of a cycle.
	GridFunction eNew{ solver.getField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

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

TEST_F(Solver1DTest, compareSpectralToBase_centered)
{
	Probes probes;
	
	maxwell::Solver solver{
	buildStandardModel(),
	probes,
	buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
	SolverOptions{}
		.setCentered()
	};
	
	maxwell::Solver solverSpectral{
	buildStandardModel(),
	probes,
	buildGaussianInitialField(E, 0.1, Vector({0.5}), unitVec(Y)),
	SolverOptions{}
		.setCentered()
		.setSpectralEO()
	};

	ASSERT_EQ(solver.getField(E, Y).Size(), solverSpectral.getField(E, Y).Size());
	ASSERT_EQ(solver.getField(H, Z).Size(), solverSpectral.getField(H, Z).Size());

	for (int i = 0; i < solver.getField(E, Y).Size(); ++i) {
		EXPECT_NEAR(solver.getField(E, Y)[i], solverSpectral.getField(E, Y)[i], 1e-3);
		EXPECT_NEAR(solver.getField(H, Z)[i], solverSpectral.getField(H, Z)[i], 1e-3);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(),1e-6);

	for (int i = 0; i < solver.getField(E, Y).Size(); ++i) {
		EXPECT_NEAR(solver.getField(E, Y)[i], solverSpectral.getField(E, Y)[i], 1e-3);
		EXPECT_NEAR(solver.getField(H, Z)[i], solverSpectral.getField(H, Z)[i], 1e-3);
	}


}

TEST_F(Solver1DTest, fieldProbeThroughSolver)
{
	Mesh m{ Mesh::MakeCartesian1D(20,5.0) };
	
	auto probes{ buildProbesWithAnExportProbe(2) };
	probes.fieldProbes = {
		FieldProbe{{2.0}}
	};

	Model model{ m, GeomTagToMaterial{}, GeomTagToBoundary{}, GeomTagToInteriorConditions{} };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.5, Vector({ 2.5 }), unitVec(Y)),
		SolverOptions{}
			.setTimeStep(5e-2)
			.setFinalTime(2.0)
			.setCentered()
	};

	solver.run();

	for (int probeNumber = 0; probeNumber < probes.fieldProbes.size(); probeNumber++) {
		std::ofstream file(getTestCaseName() + std::to_string(probeNumber) + ".txt");
		file << "Time // Ex // Ey // Ez // Hx // Hy // Hz //""\n";
		for (const auto& fm : solver.getFieldProbe(probeNumber).getFieldMovies()) {
			std::stringstream time, Ex, Ey, Ez, Hx, Hy, Hz;
			time << std::scientific << std::setprecision(7) << (fm.first);
			Ex << std::scientific << std::setprecision(7) << fm.second.Ex; Ey << std::scientific << std::setprecision(7) << fm.second.Ey; Ez << std::scientific << std::setprecision(7) << fm.second.Ez;
			Hx << std::scientific << std::setprecision(7) << fm.second.Hx; Hy << std::scientific << std::setprecision(7) << fm.second.Hy; Hz << std::scientific << std::setprecision(7) << fm.second.Hz;
			file << time.str() + " " + Ex.str() + " " + Ey.str() + " " + Ez.str() + " " + Hx.str() + " " + Hy.str() + " " + Hz.str() + "\n";
		}
	}


}

TEST_F(Solver1DTest, interior_boundary_marking_centered)
{
	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "1D_IntBdr_Line.msh"),1, 0) };
	auto probes{ buildProbesWithAnExportProbe(10) };

	GeomTagToBoundary att2Bdr{ {2,BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {3, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, att2Bdr, att2IntCond };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.30, Vector({2.0}), unitVec(Y)),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(8.0)
	};

	solver.run();

}

TEST_F(Solver1DTest, interior_boundary_marking_upwind)
{
	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "1D_IntBdr_Line.msh"),1, 0) };
	auto probes{ buildProbesWithAnExportProbe(10) };

	GeomTagToBoundary att2Bdr{ {2,BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {3, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, att2Bdr, att2IntCond };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.30, Vector({2.0}), unitVec(Y)),
		SolverOptions{}
			.setCFL(0.5)
			.setFinalTime(8.0)
	};

	solver.run();

}

TEST_F(Solver1DTest, interior_boundary_marking_centered_RtL)
{
	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "1D_IntBdr_Line.msh"),1, 0) };
	auto probes{ buildProbesWithAnExportProbe(10) };

	GeomTagToBoundary att2Bdr{ {2,BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {3, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, att2Bdr, att2IntCond};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.30, Vector({6.0}), unitVec(Y)),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(8.0)
	};

	solver.run();

}

TEST_F(Solver1DTest, interior_boundary_marking_upwind_RtL)
{
	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "1D_IntBdr_Line.msh"),1, 0) };
	auto probes{ buildProbesWithAnExportProbe(10) };

	GeomTagToBoundary att2Bdr{ {2,BdrCond::PEC} };
	GeomTagToInteriorConditions att2IntCond{ {3, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterial{}, att2Bdr, att2IntCond };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.30, Vector({6.0}), unitVec(Y)),
		SolverOptions{}
			.setCFL(0.5)
			.setFinalTime(8.0)
	};

	solver.run();

}