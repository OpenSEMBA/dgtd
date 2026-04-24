#include <gtest/gtest.h>

#include "SourceFixtures.h"

#include "solver/ProbesManager.h"
#include "solver/SourcesManager.h"
#include <solver/Solver.h>

#include <fstream>

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class ProbesManagerTest : public ::testing::Test {
	
};

TEST_F(ProbesManagerTest, exporterProbe)
{
	Mesh smesh{ Mesh::MakeCartesian1D(20, 1.0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes, fields };

	Probes ps;
	ps.exporterProbes = { ExporterProbe{"ProbesManagerTest"} };

	ProbesManager pM{ ps, fes, fields, SolverOptions{} };

	ASSERT_NO_THROW(pM.updateProbes(0.0));
}

TEST_F(ProbesManagerTest, fieldProbe)
{
	Mesh smesh{ Mesh::MakeCartesian1D(5) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes, fields };

	Probes probes;
	probes.fieldProbes = {
		FieldProbe{E, Z, {0.5}}
	};

	ProbesManager pM{ probes, fes, fields, SolverOptions{} };

	auto tol {1e-2};
	ASSERT_NO_THROW(pM.updateProbes(0.0));
	EXPECT_NO_THROW(pM.getFieldProbe(0));
	for (const auto& [t, f] : pM.getFieldProbe(0).getFieldMovie()){
		EXPECT_NEAR(1.0, f, tol);
	}

}

TEST_F(ProbesManagerTest, morStateProbe_writesTimeAsFirstLine)
{
	Mesh smesh{ Mesh::MakeCartesian1D(5, 1.0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };

	MORStateProbe morProbe;
	morProbe.name = "TestMORStateTime";
	morProbe.record_time_start = 2.5;
	morProbe.record_time_final = 2.5;
	morProbe.saves = 1;

	Probes ps;
	ps.morStateProbes = { morProbe };

	SolverOptions opts;
	opts.final_time = 2.5;

	ProbesManager pM{ ps, fes, fields, opts };
	pM.setCaseName("morStateProbeWriteTimeTest");

	// Reproduce the Solver loop pattern: initial call at t=0, then step through with dt=1e-3
	const double dt = 1e-3;
	const double final_time = opts.final_time;
	double time = 0.0;
	pM.updateProbes(time);
	while (time <= final_time - 1e-8 * dt) {
		time = std::min(time + dt, final_time);
		pM.updateProbes(time);
	}

	std::string file_path = "Exports/single-core/morStateProbeWriteTimeTest/MORStateProbes/TestMORStateTime/x_0";
	std::ifstream f(file_path);
	ASSERT_TRUE(f.is_open()) << "Expected file at: " << file_path;

	double saved_time;
	f >> saved_time;
	EXPECT_NEAR(2.5, saved_time, 1e-10);
}

TEST_F(ProbesManagerTest, morStateProbe_recordTimesConsecutive)
{
	Mesh smesh{ Mesh::MakeCartesian1D(5, 1.0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };

	MORStateProbe morProbe;
	morProbe.name = "TestMORStateConsecutive";
	morProbe.record_time_start = 0.0;
	morProbe.record_time_final = 6.0;
	morProbe.saves = 7;

	Probes ps;
	ps.morStateProbes = { morProbe };

	SolverOptions opts;
	opts.final_time = 6.0;

	ProbesManager pM{ ps, fes, fields, opts };
	pM.setCaseName("morStateProbeConsecutiveTest");

	// Reproduce the Solver loop pattern: initial call at t=0, then step through with dt=1e-3
	const double dt = 1e-3;
	const double final_time = opts.final_time;
	double time = 0.0;
	pM.updateProbes(time);
	while (time <= final_time - 1e-8 * dt) {
		time = std::min(time + dt, final_time);
		pM.updateProbes(time);
	}

	std::string base_path = "Exports/single-core/morStateProbeConsecutiveTest/MORStateProbes/TestMORStateConsecutive/";
	for (int i = 0; i < 7; ++i) {
		std::ifstream f(base_path + "x_" + std::to_string(i));
		ASSERT_TRUE(f.is_open()) << "Expected file x_" << i;

		double saved_time;
		f >> saved_time;
		EXPECT_NEAR(static_cast<double>(i), saved_time, dt)
			<< "Expected time " << i << " in file x_" << i;
	}
}
