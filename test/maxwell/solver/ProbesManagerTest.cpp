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
	Mesh mesh{ Mesh::MakeCartesian1D(20, 1.0) };
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec };
	Fields fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes };
	sM.setInitialFields(fields);

	Probes ps;
	ps.exporterProbes = { ExporterProbe{"ProbesManagerTest"} };

	ProbesManager pM{ ps, fes, fields, SolverOptions{} };

	ASSERT_NO_THROW(pM.updateProbes(0.0));
}

TEST_F(ProbesManagerTest, fieldProbe)
{
	Mesh m{ Mesh::MakeCartesian1D(5) };
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &m, &fec };
	Fields fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes };
	sM.setInitialFields(fields);

	Probes probes;
	probes.fieldProbes = {
		FieldProbe{{0.5}}
	};

	ProbesManager pM{ probes, fes, fields, SolverOptions{} };

	ASSERT_NO_THROW(pM.updateProbes(0.0));
	EXPECT_NO_THROW(pM.getFieldProbe(0));
	EXPECT_TRUE(pM.getFieldProbe(0).getFieldMovies().at(0).Ex);
}
