#include "gtest/gtest.h"
#include <math.h>

#include "maxwell/Material.h"

class TestMaxwellMaterial : public ::testing::Test {
};

TEST_F(TestMaxwellMaterial, checkImpedanceAndConductance)
{
	maxwell::Material mat1(1.0, 2.0);
	maxwell::Material mat2(100.0, 1.0);
	maxwell::Material mat3(10.0, 20.0);
	maxwell::Material mat4(9.0, 30.0);

	EXPECT_EQ(sqrt(2.0  / 1.0),   mat1.getImpedance());
	EXPECT_EQ(sqrt(1.0  / 100.0), mat2.getImpedance());
	EXPECT_EQ(sqrt(20.0 / 10.0),  mat3.getImpedance());
	EXPECT_EQ(sqrt(30.0 / 9.0),   mat4.getImpedance());
	EXPECT_EQ(sqrt(1.0  / 2.0),   mat1.getConductance());
	EXPECT_EQ(sqrt(100.0/ 1.0),   mat2.getConductance());
	EXPECT_EQ(sqrt(10.0 / 20.0),  mat3.getConductance());
	EXPECT_EQ(sqrt(9.0  / 30.0),  mat4.getConductance());

}