#include <gtest/gtest.h>
#include <cmath>

#include "math/PhysicalConstants.h"

using namespace maxwell::physicalConstants;

class PhysicalConstantsTest : public ::testing::Test {
};

TEST_F(PhysicalConstantsTest, normalizedUnitsAreUnity)
{
	EXPECT_DOUBLE_EQ(1.0, speedOfLight);
	EXPECT_DOUBLE_EQ(1.0, vacuumPermittivity);
	EXPECT_DOUBLE_EQ(1.0, vacuumPermeability);
}

TEST_F(PhysicalConstantsTest, normalizedDerivedConstants)
{
	EXPECT_DOUBLE_EQ(vacuumPermeability * speedOfLight, freeSpaceImpedance);
	EXPECT_DOUBLE_EQ(1.0 / (4.0 * M_PI * vacuumPermittivity), invFourPiEps0);
	EXPECT_DOUBLE_EQ(1.0 / (4.0 * M_PI), invFourPi);
}

TEST_F(PhysicalConstantsTest, SI_speedOfLightRelationship)
{
	// c^2 * eps0 * mu0 = 1
	double product = speedOfLight_SI * speedOfLight_SI
	               * vacuumPermittivity_SI * vacuumPermeability_SI;
	EXPECT_NEAR(1.0, product, 1e-6);
}

TEST_F(PhysicalConstantsTest, SI_freeSpaceImpedance)
{
	// Z0 = mu0 * c
	EXPECT_DOUBLE_EQ(vacuumPermeability_SI * speedOfLight_SI, freeSpaceImpedance_SI);

	// Z0 = sqrt(mu0 / eps0)
	EXPECT_NEAR(
		std::sqrt(vacuumPermeability_SI / vacuumPermittivity_SI),
		freeSpaceImpedance_SI,
		1e-3
	);
}

TEST_F(PhysicalConstantsTest, SI_invFourPiEps0)
{
	EXPECT_DOUBLE_EQ(
		1.0 / (4.0 * M_PI * vacuumPermittivity_SI),
		invFourPiEps0_SI
	);
}

TEST_F(PhysicalConstantsTest, SI_constantsArePositive)
{
	EXPECT_GT(speedOfLight_SI, 0.0);
	EXPECT_GT(vacuumPermittivity_SI, 0.0);
	EXPECT_GT(vacuumPermeability_SI, 0.0);
	EXPECT_GT(freeSpaceImpedance_SI, 0.0);
}
