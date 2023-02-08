#include "gtest/gtest.h"
#include "SourceFixtures.h"
#include "GlobalFunctions.h"

#include "maxwell/Solver.h"

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

		return Model{
			Mesh::MakeCartesian1D(numberOfElements, 1.0),
			AttributeToMaterial{},
			buildAttrToBdrMap1D(bdrL, bdrR) 
		};
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

	Probes buildProbesWithAnExportProbe()
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

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}

};

TEST_F(Solver1DTest, pec_centered)
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
		buildGaussianInitialField(E, Y),
		SolverOptions{}
			.setCentered()
	};
	
	GridFunction eOld{ solver.getFields().E[Y] };
	auto normOld{ solver.getFields().getNorml2() };
	
	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld); 
	
	solver.run();
	
	// Checks that field is almost the same as initially because the completion of a cycle.
	GridFunction eNew{ solver.getFields().E[Y] };
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

TEST_F(Solver1DTest, pmc_centered)
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
		buildStandardModel(defaultNumberOfElements, BdrCond::PMC,BdrCond::PMC),
		buildProbesWithAnExportProbe(),
		buildGaussianInitialField(H, Z),
		SolverOptions{}
			.setTimeStep(2.5e-3)
			.setCentered()
	};

	GridFunction hOld{ solver.getFields().H[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction hNew{ solver.getFields().H[Z] };

	EXPECT_NE(0.0, normOld);
	EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);
}

TEST_F(Solver1DTest, pec_upwind)
{
	maxwell::Solver solver{
		buildStandardModel(),
		buildProbesWithAnExportProbe(),
		buildGaussianInitialField(E, Y),
		SolverOptions{}
			.setCFL(0.65)
	};

	GridFunction eOld{ solver.getFields().E[Y] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction eNew{ solver.getFields().E[Y] };

	EXPECT_NE(0.0, normOld);
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-2);
}

TEST_F(Solver1DTest, pmc_upwind)
{
	maxwell::Solver solver{
		buildStandardModel(defaultNumberOfElements, BdrCond::PMC,BdrCond::PMC),
		buildProbesWithAnExportProbe(),
		buildGaussianInitialField(H, Z),
		SolverOptions{}
			.setCFL(0.65)
	};

	GridFunction hOld{ solver.getFields().H[Z] };
	auto normOld{ solver.getFields().getNorml2() };
	solver.run();
	GridFunction hNew{ solver.getFields().H[Z] };

	EXPECT_NE(0.0, normOld);
	EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), 1e-2);
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-2);
}

TEST_F(Solver1DTest, sma)
{
	maxwell::Solver solver(
		buildStandardModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
		buildProbesWithAnExportProbe(),
		buildGaussianInitialField(E, Y),
		SolverOptions{}
			.setTimeStep(2.5e-3)
			.setFinalTime(1.25)
	);

	solver.run();

	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), 2e-3);
}

TEST_F(Solver1DTest, sma_e_xyz)
{
	for (const auto& x : { X, Y, Z }) {
		Probes probes;
		probes.pointProbes = {
			PointProbe{E, x, {0.0} },
			PointProbe{E, x, {0.5} },
			PointProbe{E, x, {1.0} }
		};

		maxwell::Solver solver(
			buildStandardModel(defaultNumberOfElements, BdrCond::SMA, BdrCond::SMA),
			probes,
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

TEST_F(Solver1DTest, twoSourceWaveTwoMaterialsReflection_SMA_PEC)
{
	// Sends a wave through a material interface. 
	// and checks reflection and transmission.
	// Ref: https://en.wikipedia.org/wiki/Reflection_coefficient
	
	auto msh{ Mesh::MakeCartesian1D(100) };

	setAttributeIntervalMesh1D({ { 2, std::make_pair(0.50, 1.0) } }, msh);
	
	Material mat1{1.0, 1.0};
	Material mat2{4.0, 1.0};

	auto probes{ buildProbesWithAnExportProbe() };
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
		buildRightTravelingWaveInitialField(
			Gaussian{ 1, 0.05, 1.0, Vector({ 0.25 })} ),
		SolverOptions{}
			.setCFL(0.65)
			.setFinalTime(1.0)
	};
		
	solver.run();

	auto reflectCoeff{
		(mat2.getImpedance() - mat1.getImpedance()) /
		(mat2.getImpedance() + mat1.getImpedance())
	};
	auto transmissionCoeff{ 1 + reflectCoeff };

	auto timeTolerance{ 0.02 };
	auto fieldTolerance{ 0.01 };

	// Checks the reflected wave.
	{
		auto frame{ solver.getPointProbe(0).findFrameWithMin() };
		EXPECT_NEAR(0.75, frame.first, timeTolerance);
		EXPECT_NEAR(reflectCoeff, frame.second, fieldTolerance);
	}

	// Checks transmitted wave.
	{
		auto frame{ solver.getPointProbe(1).findFrameWithMax() };
		auto expectedTimeOfArrival{ 0.25 + 0.25 / mat2.getSpeedOfLight() };
		EXPECT_NEAR(expectedTimeOfArrival, frame.first, timeTolerance);
		EXPECT_NEAR(transmissionCoeff, frame.second, fieldTolerance);
	}

}

TEST_F(Solver1DTest, totalfield_upwind)
{
	Mesh mesh{ Mesh::LoadFromFile("./testData/verylonglineTFSF.mesh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::SMA}, {301,BdrCond::TotalFieldIn} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, AttributeToInteriorBoundary{} };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 30;

	maxwell::Solver solver{
		model,
		probes,
		buildPlaneWave(X,0.2,1.5,1.0),
		SolverOptions{}
			.setCFL(0.65)
			.setFinalTime(5.0)
			.setOrder(3)
	};

	solver.run();

	EXPECT_TRUE(false);
}

TEST_F(Solver1DTest, totalfieldin_intbdr_upwind)
{
	Mesh mesh{ Mesh::LoadFromFile("./testData/longlineIntBdr.mesh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC} };
	AttributeToInteriorBoundary attToIntBdr{ {301,BdrCond::TotalFieldIn} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 20;

	maxwell::Solver solver{
		model,
		probes,
		buildPlaneWave(X,0.2,1.5,1.0),
		SolverOptions{}
			.setCFL(0.65)
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	EXPECT_TRUE(false);
}

TEST_F(Solver1DTest, totalfieldin_intbdr_centered)
{
	Mesh mesh{ Mesh::LoadFromFile("./testData/longlineIntBdr.mesh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC} };
	AttributeToInteriorBoundary attToIntBdr{ {301,BdrCond::TotalFieldIn} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 20;

	maxwell::Solver solver{
		model,
		probes,
		buildPlaneWave(X,0.2,1.5,1.0),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	EXPECT_TRUE(false);
}

TEST_F(Solver1DTest, totalfieldinout_intbdr_centered)
{
	Mesh mesh{ Mesh::LoadFromFile("./testData/LineTFSFInOut.mesh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC} };
	AttributeToInteriorBoundary attToIntBdr{ {301,BdrCond::TotalFieldIn}, {302, BdrCond::TotalFieldOut} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 20;

	maxwell::Solver solver{
		model,
		probes,
		buildPlaneWave(X,0.2,1.5,1.0),
		SolverOptions{}
			.setCFL(0.5)
			.setCentered()
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	EXPECT_TRUE(false);
}

TEST_F(Solver1DTest, totalfieldinout_intbdr_upwind_pec)
{
	Mesh mesh{ Mesh::LoadFromFile("./testData/LineTFSFInOut.mesh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::PEC} };
	AttributeToInteriorBoundary attToIntBdr{ {301,BdrCond::TotalFieldIn}, {302, BdrCond::TotalFieldOut} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 20;

	maxwell::Solver solver{
		model,
		probes,
		buildPlaneWave(X,0.2,1.5,1.0),
		SolverOptions{}
			.setCFL(0.5)
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	EXPECT_TRUE(false);
}

TEST_F(Solver1DTest, totalfieldinout_intbdr_upwind_sma)
{
	Mesh mesh{ Mesh::LoadFromFile("./testData/LineTFSFInOut.mesh",1,0) };
	AttributeToBoundary attToBdr{ {2,BdrCond::SMA} };
	AttributeToInteriorBoundary attToIntBdr{ {301,BdrCond::TotalFieldIn}, {302, BdrCond::TotalFieldOut} };
	Model model{ mesh, AttributeToMaterial{}, attToBdr, attToIntBdr };

	auto probes{ buildProbesWithAnExportProbe() };
	probes.exporterProbes[0].visSteps = 10;

	maxwell::Solver solver{
		model,
		probes,
		buildPlaneWave(X,0.2,1.5,1.0),
		SolverOptions{}
			.setCFL(0.5)
			.setFinalTime(5.0)
			.setOrder(2)
	};

	solver.run();

	EXPECT_TRUE(false);
}
TEST_F(Solver1DTest, DISABLED_resonant_mode_upwind)
{
	// Resonant mode inside a PEC box. 
	auto probes{ buildProbesWithAnExportProbe() };

	double finalTime{ 1.2 };

	maxwell::Solver solver{
		buildStandardModel(),
		probes,
		buildResonantModeInitialField(E, Y, {1}),
		SolverOptions{}
			.setFinalTime(finalTime)
			.setCFL(0.5)
	};

	Vector eOld(solver.getFields().E[Y].Size());
	eOld = solver.getFields().E[Y];
	
	solver.run();

	GridFunction eNew{ solver.getFields().E[Y] };
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