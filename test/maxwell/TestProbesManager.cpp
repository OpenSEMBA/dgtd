#include "gtest/gtest.h"

#include "maxwell/ProbesManager.h"

using namespace maxwell;
using namespace mfem;

class TestProbesManager : public ::testing::Test {
	
};

TEST_F(TestProbesManager, exporterProbe)
{
	Mesh mesh{ Mesh::MakeCartesian1D(10, 1.0) };
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec };

	Probes ps;
	ps.exporterProbes = { ExporterProbe{"ProbesManagerTest"} };

	Vector solution{ fes.GetNDofs()*6 };

	//ProbesManager pM{ {} };
}
