#include "gtest/gtest.h"
#include <math.h>

#include "maxwell/Material.h"

using namespace maxwell;

class TestMaxwellMaterial : public ::testing::Test {
};
TEST_F(TestMaxwellMaterial, checkImpedanceAndConductance)
{
	Material mat1(1.0, 2.0);
	Material mat2(100.0, 1.0);
	Material mat3(10.0, 20.0);
	Material mat4(9.0, 30.0);

	EXPECT_EQ(sqrt(  2.0  /   1.0), mat1.getImpedance());
	EXPECT_EQ(sqrt(  1.0  / 100.0), mat2.getImpedance());
	EXPECT_EQ(sqrt( 20.0  /  10.0), mat3.getImpedance());
    EXPECT_EQ(sqrt(  1.0  /   2.0), mat1.getConductance());
	EXPECT_EQ(sqrt(100.0  /   1.0), mat2.getConductance());
	EXPECT_EQ(sqrt( 10.0  /  20.0), mat3.getConductance());


}

TEST_F(TestMaxwellMaterial, invalidEpsOrMu)
{
	EXPECT_ANY_THROW(Material( 0.5,    1.0));
	EXPECT_ANY_THROW(Material( 1.0,    0.5));
	EXPECT_ANY_THROW(Material( 20.0,  -1.0));
	EXPECT_ANY_THROW(Material(-20.0,   1.0));
	EXPECT_ANY_THROW(Material(-10.0, -10.0));
}