#include <gtest/gtest.h>
#include <math.h>

#include "components/Material.h"

using namespace maxwell;

class MaterialTest : public ::testing::Test {
};

TEST_F(MaterialTest, impedanceAndConductance)
{
	Material mat1(1.0  , 2.0 , 0.0);
	Material mat2(100.0, 1.0 , 0.0);
	Material mat3(10.0 , 20.0, 0.0);
	Material mat4(9.0  , 30.0, 0.0);

	EXPECT_EQ(sqrt(  2.0  /   1.0), mat1.getImpedance());
	EXPECT_EQ(sqrt(  1.0  / 100.0), mat2.getImpedance());
	EXPECT_EQ(sqrt( 20.0  /  10.0), mat3.getImpedance());
    EXPECT_EQ(sqrt(  1.0  /   2.0), mat1.getAdmitance());
	EXPECT_EQ(sqrt(100.0  /   1.0), mat2.getAdmitance());
	EXPECT_EQ(sqrt( 10.0  /  20.0), mat3.getAdmitance());
}

TEST_F(MaterialTest, invalidEpsOrMu)
{
	EXPECT_ANY_THROW(Material( 0.5,    1.0, 0.0));
	EXPECT_ANY_THROW(Material( 1.0,    0.5, 0.0));
	EXPECT_ANY_THROW(Material( 20.0,  -1.0, 0.0));
	EXPECT_ANY_THROW(Material(-20.0,   1.0, 0.0));
	EXPECT_ANY_THROW(Material(-10.0, -10.0, 0.0));
}

TEST_F(MaterialTest, electricalConductivity)
{
	Material vacuum{ buildVacuumMaterial() };
	EXPECT_EQ(0.0, vacuum.getConductivity());
}

TEST_F(MaterialTest, invalidSigma)
{
	EXPECT_ANY_THROW(Material(1.0, 2.0, -0.5 ));
	EXPECT_ANY_THROW(Material(2.0, 3.0, -5e-8));
}