#include "gtest/gtest.h"
#include <math.h>

#include "Material.h"

using namespace maxwell;

class MaterialTest : public ::testing::Test {
};
TEST_F(MaterialTest, impedanceAndConductance)
{
	Material mat1(1.0, 2.0);
	Material mat2(100.0, 1.0);
	Material mat3(10.0, 20.0);
	Material mat4(9.0, 30.0);

	EXPECT_EQ(sqrt(  2.0  /   1.0), mat1.getImpedance());
	EXPECT_EQ(sqrt(  1.0  / 100.0), mat2.getImpedance());
	EXPECT_EQ(sqrt( 20.0  /  10.0), mat3.getImpedance());
    EXPECT_EQ(sqrt(  1.0  /   2.0), mat1.getAdmitance());
	EXPECT_EQ(sqrt(100.0  /   1.0), mat2.getAdmitance());
	EXPECT_EQ(sqrt( 10.0  /  20.0), mat3.getAdmitance());
}

TEST_F(MaterialTest, invalidEpsOrMu)
{
	EXPECT_ANY_THROW(Material( 0.5,    1.0));
	EXPECT_ANY_THROW(Material( 1.0,    0.5));
	EXPECT_ANY_THROW(Material( 20.0,  -1.0));
	EXPECT_ANY_THROW(Material(-20.0,   1.0));
	EXPECT_ANY_THROW(Material(-10.0, -10.0));
}