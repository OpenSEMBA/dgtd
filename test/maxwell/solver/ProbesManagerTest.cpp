#include <gtest/gtest.h>

#include "SourceFixtures.h"

#include "solver/ProbesManager.h"
#include "solver/SourcesManager.h"
#include <solver/Solver.h>

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
	Fields fields{ fes };
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
	Fields fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes, fields };

	Probes probes;
	probes.pointProbes = {
		PointProbe{{0.5}}
	};

	ProbesManager pM{ probes, fes, fields, SolverOptions{} };

	ASSERT_NO_THROW(pM.updateProbes(0.0));
	EXPECT_NO_THROW(pM.getFieldProbe(0));
	EXPECT_TRUE(pM.getFieldProbe(0).getFieldMovies().at(0).Ez == 1.0);
}
