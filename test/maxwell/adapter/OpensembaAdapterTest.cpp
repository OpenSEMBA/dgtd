#include <gtest/gtest.h>

#include "adapter/OpensembaAdapter.h"
#include "TestUtils.h"

using namespace maxwell;

std::string caseFile(const std::string& caseName)
{
	return caseName + "/" + caseName + ".smb.json";
}

class OpensembaAdapterTest : public ::testing::Test {
};

TEST_F(OpensembaAdapterTest, 2D_box_resonant_mode_options)
{
	auto fn{ smbInputsFolder() + caseFile("2D_box_resonant_mode") };
	OpensembaAdapter smbAdapter(fn);

	auto opts{ smbAdapter.readSolverOptions() };

	EXPECT_EQ(0.0, opts.timeStep);
	EXPECT_EQ(2.0, opts.finalTime);
	EXPECT_EQ(0.8, opts.cfl);

	EXPECT_EQ(2, opts.evolution.order);
	EXPECT_EQ(FluxType::Centered, opts.evolution.fluxType);
	EXPECT_EQ(false, opts.evolution.spectral);
}

TEST_F(OpensembaAdapterTest, 2D_box_resonant_mode_problem)
{
	auto fn{ smbInputsFolder() + caseFile("2D_box_resonant_mode") };
	OpensembaAdapter smbAdapter(fn);

	auto prob{ smbAdapter.readProblem() };

	EXPECT_GT(0, prob.model.getMesh().GetNE());
	EXPECT_GT(0, prob.model.getMesh().GetNBE());
	EXPECT_EQ(1, prob.model.numberOfMaterials());
	EXPECT_EQ(1, prob.model.numberOfBoundaryMaterials());

	EXPECT_EQ(1, prob.probes.exporterProbes.size());
	EXPECT_EQ(0, prob.probes.pointProbes.size());

	EXPECT_EQ(1, prob.sources.size());
}