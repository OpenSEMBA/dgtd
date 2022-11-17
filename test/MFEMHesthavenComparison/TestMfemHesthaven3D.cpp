#include <gtest/gtest.h>

#include "maxwell/Types.h"
#include "TestMfemHesthavenFunctions.h"
#include "GlobalFunctions.h"


namespace maxwell{
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
	double JacobianFactor_ = 0.125;

};

TEST_F(MFEMHesthaven3D, checkNodalPositions)
{
	mesh_ = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::TETRAHEDRON);
	fec_ = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get(), 3, Ordering::byVDIM);

	GridFunction mfemNodes(fes_.get());
	mesh_.GetNodes(mfemNodes);

	mfemNodes.Print(std::cout);

	EXPECT_TRUE(false);
	
}

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

	EXPECT_TRUE(MFEMCutMass.isApprox(JacobianFactor_ * hesthavenMass, 1e-4));

}

TEST_F(MFEMHesthaven3D, checkDrOperator3D)
{
	set3DFES(2);

	Eigen::MatrixXd hesthavenDr{
		{-1.5,  2.0, -0.5,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
		{-0.5,  0.0,  0.5,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
		{ 0.5, -2.0,  1.5,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
		{-0.5,  1.0, -0.5, -1.0, 1.0, 0.0,  0.0, 0.0, 0.0, 0.0},
		{ 0.5, -1.0,  0.5, -1.0, 1.0, 0.0,  0.0, 0.0, 0.0, 0.0},
		{ 0.5,  0.0, -0.5, -2.0, 2.0, 0.0,  0.0, 0.0, 0.0, 0.0},
		{-0.5,  1.0, -0.5,  0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0},
		{ 0.5, -1.0,  0.5,  0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0},
		{ 0.5,  0.0, -0.5, -1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0},
		{ 0.5,  0.0, -0.5,  0.0, 0.0, 0.0, -2.0, 2.0, 0.0, 0.0}
	};

	auto MFEMDr{
		buildMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(Y,*fes_)
	};

	auto MFEMCutDr = MFEMDr.block<10, 10>(0, 0);

	EXPECT_TRUE(MFEMCutDr.isApprox(JacobianFactor_ * hesthavenDr, 1e-4));

}

TEST_F(MFEMHesthaven3D, DerivativeOperators_onetetra)
{
	Mesh meshManual = Mesh::LoadFromFile("./TestData/onetetra.mesh");
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	auto MFEMmass =	buildInverseMassMatrixEigen(*fesManual);
	auto MFEMSX = buildNormalStiffnessMatrixEigen(X, *fesManual);
	auto MFEMSY = buildNormalStiffnessMatrixEigen(Y, *fesManual);
	auto MFEMSZ = buildNormalStiffnessMatrixEigen(Z, *fesManual);
	auto MFEMDr = MFEMmass * MFEMSX;
	auto MFEMDs = MFEMmass * MFEMSY;
	auto MFEMDt = MFEMmass * MFEMSZ;

	Eigen::MatrixXd DrOperatorHesthaven{
	{  -0.5,  0.5, 0.0, 0.0},
	{  -0.5,  0.5, 0.0, 0.0},
	{  -0.5,  0.5, 0.0, 0.0},
	{  -0.5,  0.5, 0.0, 0.0}
	};

	Eigen::MatrixXd DsOperatorHesthaven{
	{  -0.5,  0.0, 0.5, 0.0},
	{  -0.5,  0.0, 0.5, 0.0},
	{  -0.5,  0.0, 0.5, 0.0},
	{  -0.5,  0.0, 0.5, 0.0}
	};

	Eigen::MatrixXd DtOperatorHesthaven{
	{  -0.5,  0.0, 0.0, 0.5},
	{  -0.5,  0.0, 0.0, 0.5},
	{  -0.5,  0.0, 0.0, 0.5},
	{  -0.5,  0.0, 0.0, 0.5}
	};

	std::cout << "Dr" << std::endl;
	std::cout << MFEMDr << std::endl;
	std::cout << "Ds" << std::endl;
	std::cout << MFEMDs << std::endl;
	std::cout << "Dt" << std::endl;
	std::cout << MFEMDt << std::endl;

	EXPECT_TRUE(MFEMDr.isApprox(DrOperatorHesthaven));
	EXPECT_TRUE(MFEMDs.isApprox(DsOperatorHesthaven));
	EXPECT_TRUE(MFEMDt.isApprox(DtOperatorHesthaven));
}

TEST_F(MFEMHesthaven3D, DerivativeOperators_fivetetra)
{
	Mesh meshManual = Mesh::LoadFromFile("./TestData/fivetetra.mesh");
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	auto MFEMmass = buildInverseMassMatrixEigen(*fesManual);
	auto MFEMSX = buildNormalStiffnessMatrixEigen(X, *fesManual);
	auto MFEMSY = buildNormalStiffnessMatrixEigen(Y, *fesManual);
	auto MFEMSZ = buildNormalStiffnessMatrixEigen(Z, *fesManual);
	auto MFEMDr = MFEMmass * MFEMSX;
	auto MFEMDs = MFEMmass * MFEMSY;
	auto MFEMDt = MFEMmass * MFEMSZ;

	Eigen::MatrixXd DrOperatorHesthaven{
	{  -0.5,  0.5, 0.0, 0.0},
	{  -0.5,  0.5, 0.0, 0.0},
	{  -0.5,  0.5, 0.0, 0.0},
	{  -0.5,  0.5, 0.0, 0.0}
	};

	Eigen::MatrixXd DsOperatorHesthaven{
	{  -0.5,  0.0, 0.5, 0.0},
	{  -0.5,  0.0, 0.5, 0.0},
	{  -0.5,  0.0, 0.5, 0.0},
	{  -0.5,  0.0, 0.5, 0.0}
	};

	Eigen::MatrixXd DtOperatorHesthaven{
	{  -0.5,  0.0, 0.0, 0.5},
	{  -0.5,  0.0, 0.0, 0.5},
	{  -0.5,  0.0, 0.0, 0.5},
	{  -0.5,  0.0, 0.0, 0.5}
	};

	std::cout << "Dr" << std::endl;
	std::cout << MFEMDr << std::endl;
	std::cout << "Ds" << std::endl;
	std::cout << MFEMDs << std::endl;
	std::cout << "Dt" << std::endl;
	std::cout << MFEMDt << std::endl;

	EXPECT_TRUE(MFEMDr.isApprox(DrOperatorHesthaven));
	EXPECT_TRUE(MFEMDs.isApprox(DsOperatorHesthaven));
	EXPECT_TRUE(MFEMDt.isApprox(DtOperatorHesthaven));
}
}
