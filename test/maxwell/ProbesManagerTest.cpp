#include "gtest/gtest.h"
#include "SourceFixtures.h"

#include "ProbesManager.h"
#include "SourcesManager.h"

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
