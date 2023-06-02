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

TEST_F(OpensembaAdapterTest, 2D_box_resonant_mode)
{
	auto fn{ smbInputsFolder() + caseFile("2D_box_resonant_mode") };
	OpensembaAdapter smbAdapter(fn);

	// SolverOptions parsing.
	auto opts{ smbAdapter.readSolverOptions() };

	EXPECT_EQ(0.0, opts.timeStep);
	EXPECT_EQ(2.0, opts.finalTime);
	EXPECT_EQ(0.8, opts.cfl);

	EXPECT_EQ(2, opts.evolution.order);
	EXPECT_EQ(FluxType::Centered, opts.evolution.fluxType);
	EXPECT_EQ(false, opts.evolution.spectral);

	// Problem parsing.
	auto p{ smbAdapter.readProblem() };

	EXPECT_GT(0, p.model.getMesh().GetNE());
	EXPECT_GT(0, p.model.getMesh().GetNBE());
	EXPECT_EQ(1, p.model.numberOfMaterials());
	EXPECT_EQ(1, p.model.numberOfBoundaryMaterials());

	ASSERT_EQ(1,                p.probes.exporterProbes.size());
	const ExporterProbe& ep{p.probes.exporterProbes[0]};
	EXPECT_EQ("maxwell_fields", ep.name);
	EXPECT_EQ(10,               ep.visSteps);
		
	ASSERT_EQ(1, p.sources.size());
	InitialField* src{ dynamic_cast<InitialField*>(p.sources.begin()->get()) };
	ASSERT_NE(nullptr, src);
	EXPECT_EQ(FieldType::E, src->fieldType());
	EXPECT_EQ(Source::Polarization({ 0.0, 0.0, 1.0 }), src->polarization());
	ASSERT_NE(nullptr, dynamic_cast<SinusoidalMode*>(src->magnitude()));
}