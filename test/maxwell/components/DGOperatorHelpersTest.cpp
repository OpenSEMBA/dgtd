#include <gtest/gtest.h>

#include "components/DGOperatorFactory.h"
#include "components/Sources.h"
#include "evolution/EvolutionOptions.h"
#include "math/Function.h"

using namespace maxwell;
using namespace mfem;

class DGOperatorHelpersTest : public ::testing::Test {
};

// --- altField ---

TEST_F(DGOperatorHelpersTest, altField_E_to_H)
{
	EXPECT_EQ(FieldType::H, altField(FieldType::E));
}

TEST_F(DGOperatorHelpersTest, altField_H_to_E)
{
	EXPECT_EQ(FieldType::E, altField(FieldType::H));
}

// --- Boundary coefficient maps ---

TEST_F(DGOperatorHelpersTest, PEC_coefficients)
{
	auto centPEC = bdrCentCoeff.at(BdrCond::PEC);
	EXPECT_DOUBLE_EQ(2.0, centPEC[0]);
	EXPECT_DOUBLE_EQ(0.0, centPEC[1]);

	auto upPEC = bdrUpwindCoeff.at(BdrCond::PEC);
	EXPECT_DOUBLE_EQ(2.0, upPEC[0]);
	EXPECT_DOUBLE_EQ(0.0, upPEC[1]);
}

TEST_F(DGOperatorHelpersTest, PMC_coefficients)
{
	auto centPMC = bdrCentCoeff.at(BdrCond::PMC);
	EXPECT_DOUBLE_EQ(0.0, centPMC[0]);
	EXPECT_DOUBLE_EQ(2.0, centPMC[1]);

	auto upPMC = bdrUpwindCoeff.at(BdrCond::PMC);
	EXPECT_DOUBLE_EQ(0.0, upPMC[0]);
	EXPECT_DOUBLE_EQ(2.0, upPMC[1]);
}

TEST_F(DGOperatorHelpersTest, SMA_coefficients)
{
	auto centSMA = bdrCentCoeff.at(BdrCond::SMA);
	EXPECT_DOUBLE_EQ(1.0, centSMA[0]);
	EXPECT_DOUBLE_EQ(1.0, centSMA[1]);

	auto upSMA = bdrUpwindCoeff.at(BdrCond::SMA);
	EXPECT_DOUBLE_EQ(1.0, upSMA[0]);
	EXPECT_DOUBLE_EQ(1.0, upSMA[1]);
}

TEST_F(DGOperatorHelpersTest, SGBC_coefficients)
{
	auto centSGBC = bdrCentCoeff.at(BdrCond::SGBC);
	EXPECT_DOUBLE_EQ(1.0, centSGBC[0]);
	EXPECT_DOUBLE_EQ(1.0, centSGBC[1]);

	auto upSGBC = bdrUpwindCoeff.at(BdrCond::SGBC);
	EXPECT_DOUBLE_EQ(1.0, upSGBC[0]);
	EXPECT_DOUBLE_EQ(1.0, upSGBC[1]);
}

TEST_F(DGOperatorHelpersTest, SurfaceCond_coefficients)
{
	auto centSC = bdrCentCoeff.at(BdrCond::SurfaceCond);
	EXPECT_DOUBLE_EQ(1.0, centSC[0]);
	EXPECT_DOUBLE_EQ(1.0, centSC[1]);

	auto upSC = bdrUpwindCoeff.at(BdrCond::SurfaceCond);
	EXPECT_DOUBLE_EQ(1.0, upSC[0]);
	EXPECT_DOUBLE_EQ(1.0, upSC[1]);
}

TEST_F(DGOperatorHelpersTest, sourceCoefficients)
{
	auto centTF = srcCentCoeff.at(BdrCond::TotalFieldIn);
	EXPECT_DOUBLE_EQ(1.0, centTF[0]);
	EXPECT_DOUBLE_EQ(1.0, centTF[1]);

	auto upTF = srcUpwindCoeff.at(BdrCond::TotalFieldIn);
	EXPECT_DOUBLE_EQ(1.0, upTF[0]);
	EXPECT_DOUBLE_EQ(1.0, upTF[1]);

	auto centSGBC = srcCentCoeff.at(BdrCond::SGBC);
	EXPECT_DOUBLE_EQ(1.0, centSGBC[0]);
	EXPECT_DOUBLE_EQ(1.0, centSGBC[1]);

	auto upSGBC = srcUpwindCoeff.at(BdrCond::SGBC);
	EXPECT_DOUBLE_EQ(1.0, upSGBC[0]);
	EXPECT_DOUBLE_EQ(1.0, upSGBC[1]);
}

TEST_F(DGOperatorHelpersTest, bdrCoeffCheck_centeredReturnsCorrectMap)
{
	auto centered = bdrCoeffCheck(0.0);
	EXPECT_EQ(centered.at(BdrCond::PEC), bdrCentCoeff.at(BdrCond::PEC));
	EXPECT_EQ(centered.at(BdrCond::PMC), bdrCentCoeff.at(BdrCond::PMC));
	EXPECT_EQ(centered.at(BdrCond::SMA), bdrCentCoeff.at(BdrCond::SMA));
}

TEST_F(DGOperatorHelpersTest, bdrCoeffCheck_upwindReturnsCorrectMap)
{
	auto upwind = bdrCoeffCheck(1.0);
	EXPECT_EQ(upwind.at(BdrCond::PEC), bdrUpwindCoeff.at(BdrCond::PEC));
	EXPECT_EQ(upwind.at(BdrCond::PMC), bdrUpwindCoeff.at(BdrCond::PMC));
	EXPECT_EQ(upwind.at(BdrCond::SMA), bdrUpwindCoeff.at(BdrCond::SMA));
}

TEST_F(DGOperatorHelpersTest, bdrCoeffMapsSizes)
{
	EXPECT_EQ(5u, bdrCentCoeff.size());
	EXPECT_EQ(5u, bdrUpwindCoeff.size());
	EXPECT_EQ(2u, srcCentCoeff.size());
	EXPECT_EQ(2u, srcUpwindCoeff.size());
}

// --- FieldOffsets ---

TEST_F(DGOperatorHelpersTest, FieldOffsets_Ex_local)
{
	int localSize = 10;
	int nbrSize = 5;
	FieldOffsets off(localSize, nbrSize, FieldType::E, X, true);

	// f=E(0), d=X(0), index = 3*0+0 = 0
	EXPECT_EQ(0,  off.rowStartOffset);
	EXPECT_EQ(10, off.rowEndOffset);
	EXPECT_EQ(0,  off.colStartOffset);
	EXPECT_EQ(10, off.colEndOffset); // local: only localSize cols
}

TEST_F(DGOperatorHelpersTest, FieldOffsets_Ex_nonLocal)
{
	int localSize = 10;
	int nbrSize = 5;
	FieldOffsets off(localSize, nbrSize, FieldType::E, X, false);

	EXPECT_EQ(0,  off.rowStartOffset);
	EXPECT_EQ(10, off.rowEndOffset);
	EXPECT_EQ(0,  off.colStartOffset);
	EXPECT_EQ(15, off.colEndOffset); // non-local: localSize + nbrSize
}

TEST_F(DGOperatorHelpersTest, FieldOffsets_Hz)
{
	int localSize = 8;
	int nbrSize = 4;
	// f=H(1), d=Z(2), index = 3*1+2 = 5
	FieldOffsets off(localSize, nbrSize, FieldType::H, Z, false);

	EXPECT_EQ(5 * 8,              off.rowStartOffset);
	EXPECT_EQ(5 * 8 + 8,          off.rowEndOffset);
	EXPECT_EQ(5 * (8 + 4),        off.colStartOffset);
	EXPECT_EQ(5 * (8 + 4) + 12,   off.colEndOffset);
}

TEST_F(DGOperatorHelpersTest, FieldOffsets_rowSize_matches_localBlockSize)
{
	for (int localSize = 1; localSize <= 20; localSize += 5) {
		for (auto f : {E, H}) {
			for (auto d : {X, Y, Z}) {
				FieldOffsets off(localSize, 3, f, d, true);
				EXPECT_EQ(localSize, off.rowEndOffset - off.rowStartOffset);
			}
		}
	}
}

// --- GlobalIndices ---

TEST_F(DGOperatorHelpersTest, GlobalIndices_allEntriesPopulated)
{
	GlobalIndices gi(10, 5);

	for (auto f : {E, H}) {
		for (auto d : {X, Y, Z}) {
			ASSERT_NE(nullptr, gi.offsets[f][d].get());
			int idx = 3 * f + d;
			EXPECT_EQ(idx * 10, gi.offsets[f][d]->rowStartOffset);
		}
	}
}

// --- EvolutionOptions ---

TEST_F(DGOperatorHelpersTest, EvolutionOptions_defaults)
{
	EvolutionOptions opts;

	EXPECT_EQ(EvolutionOperatorType::Global, opts.op);
	EXPECT_EQ(2, opts.order);
	EXPECT_DOUBLE_EQ(1.0, opts.alpha);
	EXPECT_FALSE(opts.spectral);
	EXPECT_FALSE(opts.export_evolution_operator);
	EXPECT_FALSE(opts.is_sgbc_solver);
}

// --- SGBCProperties ---

TEST_F(DGOperatorHelpersTest, SGBCProperties_empty)
{
	SGBCProperties props;

	EXPECT_DOUBLE_EQ(0.0, props.totalWidth());
	EXPECT_EQ(0u, props.totalSegments());
	EXPECT_EQ(1u, props.maxOrder());
}

TEST_F(DGOperatorHelpersTest, SGBCProperties_singleLayer)
{
	SGBCProperties props;
	props.layers.emplace_back(Material(2.0, 1.0, 0.0), 0.5);
	props.layers.back().num_of_segments = 20;
	props.layers.back().order = 3;

	EXPECT_DOUBLE_EQ(0.5, props.totalWidth());
	EXPECT_EQ(20u, props.totalSegments());
	EXPECT_EQ(3u,  props.maxOrder());
}

TEST_F(DGOperatorHelpersTest, SGBCProperties_multipleLayers)
{
	SGBCProperties props;
	props.layers.emplace_back(Material(2.0, 1.0, 0.0), 0.3);
	props.layers.back().num_of_segments = 10;
	props.layers.back().order = 2;

	props.layers.emplace_back(Material(5.0, 1.0, 50.0), 0.7);
	props.layers.back().num_of_segments = 30;
	props.layers.back().order = 4;

	EXPECT_DOUBLE_EQ(1.0, props.totalWidth());
	EXPECT_EQ(40u, props.totalSegments());
	EXPECT_EQ(4u,  props.maxOrder());
}

TEST_F(DGOperatorHelpersTest, SGBCLayer_defaultValues)
{
	SGBCLayer layer(Material(1.0, 1.0, 0.0), 0.1);

	EXPECT_EQ(10u, layer.num_of_segments);
	EXPECT_EQ(1u,  layer.order);
	EXPECT_DOUBLE_EQ(0.0, layer.n_skin_depths);
	EXPECT_DOUBLE_EQ(0.1, layer.width);
}

// --- Sources container ---

TEST_F(DGOperatorHelpersTest, SourcesContainer_emptyByDefault)
{
	Sources srcs;
	EXPECT_EQ(0u, srcs.size());
}

TEST_F(DGOperatorHelpersTest, SourcesContainer_addAndSize)
{
	Sources srcs;
	mfem::Vector mean({0.0});
	Gaussian g(1.0, mean, 1);
	mfem::Vector pol({1.0, 0.0, 0.0});
	mfem::Vector center({0.0, 0.0, 0.0});

	srcs.add(std::make_unique<InitialField>(g, FieldType::E, pol, center));
	EXPECT_EQ(1u, srcs.size());

	srcs.add(std::make_unique<InitialField>(g, FieldType::H, pol, center));
	EXPECT_EQ(2u, srcs.size());
}

TEST_F(DGOperatorHelpersTest, SourcesContainer_addByReference)
{
	Sources srcs;
	mfem::Vector mean({0.0});
	Gaussian g(1.0, mean, 1);
	mfem::Vector pol({1.0, 0.0, 0.0});
	mfem::Vector center({0.0, 0.0, 0.0});

	InitialField src(g, FieldType::E, pol, center);
	auto* ptr = srcs.add(src);
	EXPECT_NE(nullptr, ptr);
	EXPECT_EQ(1u, srcs.size());
}

TEST_F(DGOperatorHelpersTest, SourcesContainer_copyConstruction)
{
	Sources srcs;
	mfem::Vector mean({0.0});
	Gaussian g(1.0, mean, 1);
	mfem::Vector pol({1.0, 0.0, 0.0});
	mfem::Vector center({0.0, 0.0, 0.0});
	srcs.add(std::make_unique<InitialField>(g, FieldType::E, pol, center));

	Sources copy(srcs);
	EXPECT_EQ(1u, copy.size());
	EXPECT_EQ(1u, srcs.size()); // original unaffected
}

TEST_F(DGOperatorHelpersTest, SourcesContainer_iteration)
{
	Sources srcs;
	mfem::Vector mean({0.0});
	Gaussian g(1.0, mean, 1);
	mfem::Vector pol({1.0, 0.0, 0.0});
	mfem::Vector center({0.0, 0.0, 0.0});
	srcs.add(std::make_unique<InitialField>(g, FieldType::E, pol, center));
	srcs.add(std::make_unique<InitialField>(g, FieldType::H, pol, center));

	int count = 0;
	for (const auto& s : srcs) {
		EXPECT_NE(nullptr, s.get());
		count++;
	}
	EXPECT_EQ(2, count);
}
