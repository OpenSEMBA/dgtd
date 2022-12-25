#include "gtest/gtest.h"

#include "mfem.hpp"

using namespace mfem;

TEST(MFEM, IntegrationRule_tetrahedron)
{
	IntegrationRule ir{ IntRules.Get(Geometry::TETRAHEDRON, 4) };

	// Integrates over constant volume. Rule is defined over a reference tetrahedron.
	auto integral{ 0.0 };
	for (auto p{ 0 }; p < ir.GetNPoints(); p++) {
		integral += ir.IntPoint(p).weight;
	}

	// Compares with volume of reference tetrahedron.
	EXPECT_EQ(1.0 / 6.0, integral);
}