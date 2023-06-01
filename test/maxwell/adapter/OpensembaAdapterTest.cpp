#include <gtest/gtest.h>

#include "adapter/OpensembaAdapter.h"
#include "TestUtils.h"

using namespace maxwell;

class OpensembaAdapterTest : public ::testing::Test {
};

TEST_F(OpensembaAdapterTest, 2D_box_resonant_mode)
{
	auto fn{ smbInputsFolder() + getCaseName() + ".smb.json" };
	OpensembaAdapter smbAdapter(fn);

	auto prob{ smbAdapter.readProblem() };
	auto opts{ smbAdapter.readSolverOptions() };

	EXPECT_GT(0, prob.model.getMesh().GetNE());
	EXPECT_GT(0, prob.model.getMesh().GetNBE());
	EXPECT_EQ(1, prob.model.numberOfMaterials());
	EXPECT_EQ(1, prob.model.numberOfBoundaryMaterials());

	EXPECT_EQ(1, prob.probes.exporterProbes.size());
	EXPECT_EQ(0, prob.probes.pointProbes.size());

	EXPECT_EQ(1, prob.sources.size());
}