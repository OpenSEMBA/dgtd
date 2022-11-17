#include <gtest/gtest.h>

#include "TestMfemHesthavenFunctions.h"
#include "GlobalFunctions.h"

using namespace mfem;

class MFEMHesthaven3D : public ::testing::Test {
protected:

	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::TETRAHEDRON);
		fec_ = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void set3DFES(
		const int order, 
		const int xElem = 1,
		const int yElem = 1, 
		const int zElem = 1,
		Element::Type eType = Element::Type::TETRAHEDRON)
	{
		mesh_ = Mesh::MakeCartesian3D(xElem, yElem, zElem, eType);
		fec_ = std::make_unique<DG_FECollection>(order, 3, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};

TEST_F(MFEMHesthaven3D, checkMassOperator3D)
{
	set3DFES(2);

	Eigen::MatrixXd hesthavenMass{
	{ 0.019048, -0.012698, 0.0031746, -0.012698, -0.019048, 0.0031746, -0.012698, -0.019048, -0.019048, 0.0031746},
    {-0.012698,	  0.10159, -0.012698,  0.050794,  0.050794, -0.019048,  0.050794,  0.050794,  0.025397, -0.019048},
    {0.0031746, -0.012698,  0.019048, -0.019048, -0.012698, 0.0031746, -0.019048, -0.012698, -0.019048, 0.0031746},
    {-0.012698,	 0.050794, -0.019048,   0.10159,  0.050794, -0.012698,  0.050794,  0.025397,  0.050794, -0.019048},
    {-0.019048,	 0.050794, -0.012698,  0.050794,   0.10159, -0.012698,  0.025397,  0.050794,  0.050794, -0.019048},
    {0.0031746, -0.019048, 0.0031746, -0.012698, -0.012698,  0.019048, -0.019048, -0.019048, -0.012698, 0.0031746},
    {-0.012698,	 0.050794, -0.019048,  0.050794,  0.025397, -0.019048,   0.10159,  0.050794,  0.050794, -0.012698},
    {-0.019048,	 0.050794, -0.012698,  0.025397,  0.050794, -0.019048,  0.050794,   0.10159,  0.050794, -0.012698},
    {-0.019048,	 0.025397, -0.019048,  0.050794,  0.050794, -0.012698,  0.050794,  0.050794,   0.10159, -0.012698},
    {0.0031746 ,-0.019048, 0.0031746, -0.019048, -0.019048, 0.0031746, -0.012698, -0.012698, -0.012698,  0.019048}
	};

	auto MFEMmass = buildMassMatrixEigen(*fes_);
	auto MFEMCutMass = MFEMmass.block<10, 10>(0, 0);

	EXPECT_TRUE(MFEMCutMass.isApprox(0.125 * hesthavenMass, 1e-4));


	
}
