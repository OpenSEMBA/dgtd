#include "components/Spherical.h"
#include "gtest/gtest.h"

using namespace maxwell;

class SphericalTest : public ::testing::Test {
};

TEST_F(SphericalTest, initSphericalVector) 
{
	std::vector<double> vector_2D({ 1.0, 2.0, 0.0 });
	std::vector<double> vector_3D({ 1.0, 1.0, 1.0 });
	mfem::Vector vector_3D_mfem_neg({ -1.0, -1.0, -1.0 });
	std::vector<double> vector_3D_bigger({ 4.4, -2.3, 9.2 });

	SphericalVector sph_2D(vector_2D);
	SphericalVector sph_3D(vector_3D);
	SphericalVector sph_3D_mfem_neg(vector_3D_mfem_neg);
	SphericalVector sph_3D_bigger(vector_3D_bigger);

	SphericalVector exp_sph_2D(2.2361, 1.5708, 1.1071);
	SphericalVector exp_sph_3D(1.7321, 0.9553, 0.7854);
	SphericalVector exp_sph_3D_mfem_neg(1.7321, 2.1863, -2.356);
	SphericalVector exp_sph_3D_bigger(10.454, 0.4949, -0.4817);

	auto tol{ 1e-3 };
	EXPECT_NEAR(exp_sph_2D.radius, sph_2D.radius, tol);
	EXPECT_NEAR(exp_sph_2D.theta, sph_2D.theta, tol);
	EXPECT_NEAR(exp_sph_2D.phi, sph_2D.phi, tol);

	EXPECT_NEAR(exp_sph_3D.radius, sph_3D.radius, tol);
	EXPECT_NEAR(exp_sph_3D.theta, sph_3D.theta, tol);
	EXPECT_NEAR(exp_sph_3D.phi, sph_3D.phi, tol);

	EXPECT_NEAR(exp_sph_3D_mfem_neg.radius, sph_3D_mfem_neg.radius, tol);
	EXPECT_NEAR(exp_sph_3D_mfem_neg.theta, sph_3D_mfem_neg.theta, tol);
	EXPECT_NEAR(exp_sph_3D_mfem_neg.phi, sph_3D_mfem_neg.phi, tol);

	EXPECT_NEAR(exp_sph_3D_bigger.radius, sph_3D_bigger.radius, tol);
	EXPECT_NEAR(exp_sph_3D_bigger.theta, sph_3D_bigger.theta, tol);
	EXPECT_NEAR(exp_sph_3D_bigger.phi, sph_3D_bigger.phi, tol);

}

TEST_F(SphericalTest, convertFieldToCartesian)
{
	std::vector<double> vector_3D({ 1.0, 1.0, 1.0 });
	SphericalVector vec_3D(vector_3D);
	auto tol{ 1e-3 };

	double Ar = 2.0;
	double At = 3.0;
	double Ap = 0.0;

	double expected_EX = 2.379445409770841;
	double expected_EY = 2.379445409770841;
	double expected_EZ = -1.2947892044039262;

	auto res_E = vec_3D.convertSphericalVectorFieldToCartesian(Ar, At, Ap);
	
	EXPECT_NEAR(expected_EX, res_E[0], tol);
	EXPECT_NEAR(expected_EY, res_E[1], tol);
	EXPECT_NEAR(expected_EZ, res_E[2], tol);

	Ar = 0.0;
	At = 0.0;
	Ap = 4.0;

	double expected_HX = -2.8284271247461903;
	double expected_HY = 2.8284271247461903;
	double expected_HZ = 0.0;

	auto res_H = vec_3D.convertSphericalVectorFieldToCartesian(Ar, At, Ap);

	EXPECT_NEAR(expected_HX, res_H[0], tol);
	EXPECT_NEAR(expected_HY, res_H[1], tol);
	EXPECT_NEAR(expected_HZ, res_H[2], tol);

	vector_3D = { -0.3, 0.5, -0.726544 };
	vec_3D = SphericalVector(vector_3D);

	Ar = 2.0;
	At = 3.0;
	Ap = 0.0;

	expected_EX = 0.5596985080849165;
	expected_EY = -0.9328308468081946;
	expected_EZ = -3.43752297320187;

	res_E = vec_3D.convertSphericalVectorFieldToCartesian(Ar, At, Ap);

	EXPECT_NEAR(expected_EX, res_E[0], tol);
	EXPECT_NEAR(expected_EY, res_E[1], tol);
	EXPECT_NEAR(expected_EZ, res_E[2], tol);

	Ar = 0.0;
	At = 0.0;
	Ap = 4.0;

	expected_HX = -3.4299717028501773;
	expected_HY = -2.0579830217101054;
	expected_HZ = 0.0;

	res_H = vec_3D.convertSphericalVectorFieldToCartesian(Ar, At, Ap);

	EXPECT_NEAR(expected_HX, res_H[0], tol);
	EXPECT_NEAR(expected_HY, res_H[1], tol);
	EXPECT_NEAR(expected_HZ, res_H[2], tol);
}