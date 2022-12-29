#include <gtest/gtest.h>

#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/Types.h"
#include "TestMfemHesthavenFunctions.h"
#include "GlobalFunctions.h"


using namespace maxwell;
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

	Eigen::Matrix<double, 4, 4> rotatorO1_{
		{0,0,0,1},
		{1,0,0,0},
		{0,1,0,0},
		{0,0,1,0}
	};

};

TEST_F(MFEMHesthaven3D, DISABLED_checkNodalPositions)
{
	mesh_ = Mesh::LoadFromFile("./testData/onetetra.mesh");
	fec_ = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get(), 3, Ordering::byVDIM);

	GridFunction mfemNodes(fes_.get());
	mesh_.GetNodes(mfemNodes);

	std::cout << mfemNodes << std::endl;

	Eigen::MatrixXd mfemEigenNodes(4, 3);
	for (Eigen::Index i = 0; i < mfemEigenNodes.rows(); i++) {
		for (Eigen::Index j = 0; j < mfemEigenNodes.cols(); j++) {
			mfemEigenNodes(i, j) = mfemNodes.Elem((i * mfemEigenNodes.cols()) + j);
		}
	}

	auto mfemEigenRotatedNodes = rotatorO1_ * mfemEigenNodes;

	Eigen::MatrixXd expectedNodes{
		{0,0,0},
		{1,0,0},
		{0,1,0},
		{0,0,1}
	};

	std::cout << rotatorO1_.transpose() * expectedNodes  << std::endl;

	EXPECT_TRUE(expectedNodes.isApprox(mfemEigenRotatedNodes, tol_));
	
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

TEST_F(MFEMHesthaven3D, DISABLED_checkDrOperator3D)
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
	auto MFEMDr = MFEMmass * MFEMSZ;
	auto MFEMDs = MFEMmass * MFEMSY;
	auto MFEMDt = MFEMmass * MFEMSX;

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

	auto DrRotated = rotatorO1_.transpose() * DrOperatorHesthaven * rotatorO1_;
	auto DsRotated = rotatorO1_.transpose() * DsOperatorHesthaven * rotatorO1_;
	auto DtRotated = rotatorO1_.transpose() * DtOperatorHesthaven * rotatorO1_;

	std::cout << "Dr" << std::endl;
	std::cout << MFEMDr << std::endl;
	std::cout << DrRotated << std::endl;
	std::cout << "Ds" << std::endl;
	std::cout << MFEMDs << std::endl;
	std::cout << DsRotated << std::endl;
	std::cout << "Dt" << std::endl;
	std::cout << MFEMDt << std::endl;
	std::cout << DtRotated << std::endl;

	EXPECT_TRUE(MFEMDr.isApprox(DrRotated));
	EXPECT_TRUE(MFEMDs.isApprox(DsRotated));
	EXPECT_TRUE(MFEMDt.isApprox(DtRotated));
}

TEST_F(MFEMHesthaven3D, DISABLED_DerivativeOperators_fivetetra)
{
	Mesh meshManual = Mesh::LoadFromFile("./testData/fivetetra.mesh");
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	auto MFEMmass = buildInverseMassMatrixEigen(*fesManual);
	auto MFEMSX = buildNormalStiffnessMatrixEigen(Z, *fesManual);
	auto MFEMSY = buildNormalStiffnessMatrixEigen(Y, *fesManual);
	auto MFEMSZ = buildNormalStiffnessMatrixEigen(X, *fesManual);
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

TEST_F(MFEMHesthaven3D, faceChecker)
{

	Mesh meshManual = Mesh::LoadFromFile("./TestData/onetetra.mesh");
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	{
		BoundaryMarker bdrMarker{ meshManual.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[0] = 1;

		auto form = std::make_unique<BilinearForm>(fesManual.get());
		form->AddBdrFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(std::vector<Direction>{X}, 1.0), bdrMarker);
		form->Assemble();
		form->Finalize();

		std::cout << "Face 0" << std::endl;
		std::cout << toEigen(*form.get()->SpMat().ToDenseMatrix()) << std::endl;
	}

	{
		BoundaryMarker bdrMarker{ meshManual.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[1] = 1;

		auto form = std::make_unique<BilinearForm>(fesManual.get());
		form->AddBdrFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(std::vector<Direction>{X}, 1.0), bdrMarker);
		form->Assemble();
		form->Finalize();

		std::cout << "Face 1" << std::endl;
		std::cout << toEigen(*form.get()->SpMat().ToDenseMatrix()) << std::endl;
	}

	{
		BoundaryMarker bdrMarker{ meshManual.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[2] = 1;

		auto form = std::make_unique<BilinearForm>(fesManual.get());
		form->AddBdrFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(std::vector<Direction>{X}, 1.0), bdrMarker);
		form->Assemble();
		form->Finalize();

		std::cout << "Face 2" << std::endl;
		std::cout << toEigen(*form.get()->SpMat().ToDenseMatrix()) << std::endl;
	}

	{
		BoundaryMarker bdrMarker{ meshManual.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[3] = 1;

		auto form = std::make_unique<BilinearForm>(fesManual.get());
		form->AddBdrFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(std::vector<Direction>{X}, 1.0), bdrMarker);
		form->Assemble();
		form->Finalize();

		std::cout << "Face 3" << std::endl;
		std::cout << toEigen(*form.get()->SpMat().ToDenseMatrix()) << std::endl;
	}

	{
		Mesh meshTriang = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE);
		std::unique_ptr<FiniteElementCollection> fecTriang = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
		std::unique_ptr<FiniteElementSpace> fesTriang = std::make_unique<FiniteElementSpace>(&meshTriang, fecTriang.get());
		auto form = std::make_unique<BilinearForm>(fesTriang.get());
		auto one = ConstantCoefficient(1.0);
		form->AddDomainIntegrator(new MassIntegrator(one));
		form->Assemble();
		form->Finalize();
		std::cout << "Mass Triang" << std::endl;
		std::cout << toEigen(*form.get()->SpMat().ToDenseMatrix()) << std::endl;
	}


}

