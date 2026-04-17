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

TEST_F(MaterialTest, speedOfWave)
{
	Material vac(1.0, 1.0, 0.0);
	EXPECT_DOUBLE_EQ(1.0, vac.getSpeedOfWave()); // 1/sqrt(1*1)

	Material mat1(4.0, 1.0, 0.0);
	EXPECT_DOUBLE_EQ(0.5, mat1.getSpeedOfWave()); // 1/sqrt(4*1)

	Material mat2(1.0, 4.0, 0.0);
	EXPECT_DOUBLE_EQ(0.5, mat2.getSpeedOfWave()); // 1/sqrt(1*4)

	Material mat3(2.0, 8.0, 0.0);
	EXPECT_DOUBLE_EQ(1.0 / sqrt(16.0), mat3.getSpeedOfWave()); // 0.25
}

TEST_F(MaterialTest, conductiveMaterialThrows)
{
	Material cond(1.0, 1.0, 100.0);
	EXPECT_EQ(100.0, cond.getConductivity());
	EXPECT_ANY_THROW(cond.getImpedance());
	EXPECT_ANY_THROW(cond.getAdmitance());
	EXPECT_ANY_THROW(cond.getSpeedOfWave());
}

TEST_F(MaterialTest, buildVacuumMaterialProperties)
{
	Material vac = buildVacuumMaterial();
	EXPECT_DOUBLE_EQ(1.0, vac.getPermittivity());
	EXPECT_DOUBLE_EQ(1.0, vac.getPermeability());
	EXPECT_DOUBLE_EQ(0.0, vac.getConductivity());
	EXPECT_DOUBLE_EQ(1.0, vac.getImpedance());
	EXPECT_DOUBLE_EQ(1.0, vac.getAdmitance());
	EXPECT_DOUBLE_EQ(1.0, vac.getSpeedOfWave());
}

TEST_F(MaterialTest, impedanceTimesAdmittanceIsOne)
{
	Material mat1(3.0,   7.0, 0.0);
	Material mat2(100.0, 1.0, 0.0);
	Material mat3(1.0,  50.0, 0.0);

	EXPECT_NEAR(1.0, mat1.getImpedance() * mat1.getAdmitance(), 1e-15);
	EXPECT_NEAR(1.0, mat2.getImpedance() * mat2.getAdmitance(), 1e-15);
	EXPECT_NEAR(1.0, mat3.getImpedance() * mat3.getAdmitance(), 1e-15);
}

TEST_F(MaterialTest, speedOfWaveRelationToImpedance)
{
	// For lossless: v = 1/sqrt(mu*eps), Z = sqrt(mu/eps)
	// So v * Z = 1/eps, v / Z = 1/mu
	Material mat(4.0, 9.0, 0.0);
	double v = mat.getSpeedOfWave();
	double Z = mat.getImpedance();

	EXPECT_NEAR(1.0 / mat.getPermittivity(), v * Z, 1e-15);
	EXPECT_NEAR(1.0 / mat.getPermeability(), v / Z, 1e-15);
}