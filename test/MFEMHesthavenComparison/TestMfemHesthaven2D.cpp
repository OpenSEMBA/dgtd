#include <gtest/gtest.h>

#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/Types.h"
#include "TestMfemHesthavenFunctions.h"
#include "GlobalFunctions.h"
#include "maxwell/MaxwellDefs.h"
#include "maxwell/MaxwellDefs1D.h"


using namespace mfem;
using namespace maxwell;

class MFEMHesthaven2D : public ::testing::Test {
protected:

	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE);
		fec_ = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES(const int order, const int nx = 1, const int ny = 1)
	{
		mesh_ = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE);
		fec_ = std::make_unique<DG_FECollection>(order, 2, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};

TEST_F(MFEMHesthaven2D, massMatrix2D)
{
	Eigen::MatrixXd hesthavenMass{
		{0.33333, 0.16667, 0.16667,     0.0,     0.0,     0.0},
		{0.16667, 0.33333, 0.16667,     0.0,     0.0,     0.0},
		{0.16667, 0.16667, 0.33333,     0.0,     0.0,     0.0},
		{    0.0,     0.0,     0.0, 0.33333, 0.16667, 0.16667},
		{    0.0,     0.0,     0.0, 0.16667, 0.33333, 0.16667},
		{    0.0,     0.0,     0.0, 0.16667, 0.16667, 0.33333}
	};

	const double scaleFactor = 0.25;

	auto MFEMMass = buildMassMatrixEigen(*fes_);

	std::cout << MFEMMass << std::endl;

	EXPECT_TRUE(MFEMMass.isApprox(hesthavenMass * scaleFactor,1e-2));
}

TEST_F(MFEMHesthaven2D, DOperators2D)
{
	std::cout << buildMassMatrixEigen(*fes_) << std::endl;
	std::cout << buildInverseMassMatrixEigen(*fes_) << std::endl;
	std::cout << buildNormalStiffnessMatrixEigen(X, *fes_) << std::endl;
	std::cout << buildNormalStiffnessMatrixEigen(Y, *fes_) << std::endl;

	Eigen::MatrixXd DrOperatorMFEM{
		buildInverseMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(X,*fes_)
	};

	Eigen::MatrixXd DrOperatorHesthaven{
		{ -0.5, 0.5, 0.0,  0.0, 0.0, 0.0},
		{ -0.5, 0.5, 0.0,  0.0, 0.0, 0.0},
		{ -0.5, 0.5, 0.0,  0.0, 0.0, 0.0},
		{  0.0, 0.0, 0.0, -0.5, 0.5, 0.0},
		{  0.0, 0.0, 0.0, -0.5, 0.5, 0.0},
		{  0.0, 0.0, 0.0, -0.5, 0.5, 0.0}
	};

	Eigen::MatrixXd DsOperatorMFEM{
		buildInverseMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(Y,*fes_)
	};

	Eigen::MatrixXd DsOperatorHesthaven{
		{ -0.5, 0.0, 0.5,  0.0, 0.0, 0.0},
		{ -0.5, 0.0, 0.5,  0.0, 0.0, 0.0},
		{ -0.5, 0.0, 0.5,  0.0, 0.0, 0.0},
		{  0.0, 0.0, 0.0, -0.5, 0.0, 0.5},
		{  0.0, 0.0, 0.0, -0.5, 0.0, 0.5},
		{  0.0, 0.0, 0.0, -0.5, 0.0, 0.5}
	};

	const double scaleFactor = 0.25;

	std::cout << DrOperatorMFEM << std::endl;
	std::cout << DsOperatorMFEM << std::endl;

	EXPECT_TRUE(DrOperatorMFEM.isApprox(scaleFactor * DrOperatorHesthaven, tol_));
	EXPECT_TRUE(DsOperatorMFEM.isApprox(scaleFactor * DsOperatorHesthaven, tol_));

}

TEST_F(MFEMHesthaven2D, manualMeshComparison)
{
	Mesh meshManual = Mesh::LoadFromFile("./TestData/twotriang.mesh", 1, 1);
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 2);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	Mesh meshAuto = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE, true);
	std::unique_ptr<FiniteElementCollection> fecAuto = std::make_unique<DG_FECollection>(1, 2);
	std::unique_ptr<FiniteElementSpace> fesAuto = std::make_unique<FiniteElementSpace>(&meshAuto, fecAuto.get());

	ASSERT_TRUE(buildMassMatrixEigen(*fesManual).isApprox(buildMassMatrixEigen(*fesAuto),tol_));
	ASSERT_TRUE(buildNormalStiffnessMatrixEigen(X, *fesManual).isApprox(buildNormalStiffnessMatrixEigen(X, *fesAuto), tol_));
	ASSERT_TRUE(buildNormalStiffnessMatrixEigen(Y, *fesManual).isApprox(buildNormalStiffnessMatrixEigen(Y, *fesAuto), tol_));

	
}
