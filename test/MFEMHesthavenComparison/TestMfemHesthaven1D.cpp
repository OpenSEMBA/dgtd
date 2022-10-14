#include <gtest/gtest.h>

#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/Types.h"
#include "TestMfemHesthavenFunctions.h"
#include "GlobalFunctions.h"
#include "maxwell/MaxwellDefs.h"
#include "maxwell/MaxwellDefs1D.h"


using namespace mfem;
using namespace maxwell;

class MFEMHesthaven1D : public ::testing::Test {
protected:

	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian1D(1);
		fec_ = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES(const int order, const int elements = 1)
	{
		mesh_ = Mesh::MakeCartesian1D(elements);
		fec_ = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};

TEST_F(MFEMHesthaven1D, MassMatrix_O1)
{
	Eigen::MatrixXd expected{
		{2.0 / 6.0, 1.0 / 6.0},
		{1.0 / 6.0, 2.0 / 6.0} 
	};

	EXPECT_TRUE(buildMassMatrixEigen(*fes_).isApprox(expected, tol_));
}

TEST_F(MFEMHesthaven1D, InverseMassMatrix_O1)
{	
	Eigen::MatrixXd expected{
		{ 2.0, -1.0},
		{-1.0,  2.0}
	};

	EXPECT_TRUE(buildInverseMassMatrixEigen(*fes_).isApprox(expected, tol_));
}

TEST_F(MFEMHesthaven1D, StiffnessMatrix_O1)
{
	Eigen::MatrixXd expected{
			{-0.5, 0.5},
			{-0.5, 0.5} 
	};

	EXPECT_TRUE(buildStiffnessMatrixEigen(*fes_).isApprox(expected, tol_));
}

TEST_F(MFEMHesthaven1D, DOperator_O1)
{
	Eigen::MatrixXd D{ buildInverseMassMatrixEigen(*fes_) * buildStiffnessMatrixEigen(*fes_) };

	Eigen::MatrixXd expected{
		{-0.5, 0.5},
		{-0.5, 0.5} 
	};

	EXPECT_TRUE(D.isApprox(expected, tol_));
}

TEST_F(MFEMHesthaven1D, DOperator_O2)
{
	setFES(2);
	Eigen::MatrixXd D{ 
		buildInverseMassMatrixEigen(*fes_) * buildStiffnessMatrixEigen(*fes_) 
	};

	Eigen::MatrixXd expected{
		{-1.5, 2.0,-0.5},
		{-0.5, 0.0, 0.5},
		{ 0.5,-2.0, 1.5} 
	};

	EXPECT_TRUE(D.isApprox(expected, tol_));
}

TEST_F(MFEMHesthaven1D, DOperator_O4)
{
	setFES(4);
	Eigen::MatrixXd D{
		buildInverseMassMatrixEigen(*fes_) * buildStiffnessMatrixEigen(*fes_)
	};

	Eigen::MatrixXd expected{
		{-5.000000000000000e+00, 6.756502488724238e+00,-2.666666666666666e+00, 1.410164177942427e+00,-5.000000000000000e-01},
		{-1.240990253030983e+00, 0.000000000000000e+00, 1.745743121887938e+00,-7.637626158259730e-01, 2.590097469690174e-01},
		{ 3.750000000000002e-01,-1.336584577695454e+00, 0.000000000000000e+00, 1.336584577695453e+00,-3.750000000000002e-01},
		{-2.590097469690172e-01, 7.637626158259738e-01,-1.745743121887938e+00, 0.000000000000000e+00, 1.240990253030984e+00},
		{ 5.000000000000000e-01,-1.410164177942426e+00, 2.666666666666663e+00,-6.756502488724239e+00, 5.000000000000001e+00} 
	};

	EXPECT_TRUE(D.isApprox(expected,tol_));
}

TEST_F(MFEMHesthaven1D, MFOperator)
{
	setFES(2, 4);
	Eigen::MatrixXd MFField_MFEM_Ones{
		buildInverseMassMatrixEigen(*fes_) * buildNormalSMAFluxOperator1D(*fes_, std::vector<int>{0}) * Eigen::VectorXd::Ones(buildInverseMassMatrixEigen(*fes_).cols())
	};

	Eigen::VectorXd ManualFieldVector{ {5.69e-05,0.0045202,0.10288,0.10288,0.67085,1.2533,1.2533,0.67085,0.10288,0.10288,0.0045202,5.69e-05} };
	Eigen::MatrixXd MFField_MFEM_Manual{
		buildInverseMassMatrixEigen(*fes_)* buildNormalSMAFluxOperator1D(*fes_, std::vector<int>{0}) * ManualFieldVector
	};

	Eigen::VectorXd MFField_Hesthaven{{-18.0, 3.0, -6.0, 0, 0, 0, 0, 0, 0, -6.0, 3.0, -18.0}};
	Eigen::VectorXd MFField_Hesthaven_Manual{ {-0.0010242,0.0001707,-0.0003414,0,0,0,0,0,0,-0.0003414,0.0001707,-0.0010242} };

	EXPECT_TRUE(MFField_MFEM_Ones.isApprox(MFField_Hesthaven,tol_));
	EXPECT_TRUE(MFField_MFEM_Manual.isApprox(MFField_Hesthaven_Manual, tol_));
} 

TEST_F(MFEMHesthaven1D, MSOperator)
{
	setFES(2, 4);
	auto res = toEigen(*buildByMult(
		*buildInverseMassMatrix(E, Model(mesh_, AttributeToMaterial{}, { {1, BdrCond::SMA}, {2, BdrCond::SMA} }), *fes_), 
		*buildDerivativeOperator(X, *fes_), *fes_)
		.get()->SpMat().ToDenseMatrix());
	auto expected = buildMatrixForMSTest();
	EXPECT_TRUE(res.isApprox(expected,tol_));
}

TEST_F(MFEMHesthaven1D, MPOperator)
{
	setFES(2, 4);
	Eigen::MatrixXd MPField_MFEM{
		buildInverseMassMatrixEigen(*fes_) * buildSMAPenaltyOperator1D(*fes_) * Eigen::VectorXd::Ones(buildInverseMassMatrixEigen(*fes_).cols())
	};

	Eigen::VectorXd MPField_Hesthaven{ {-18.0, 3.0, -6.0, 0, 0, 0, 0, 0, 0, -6.0, 3.0, -18.0} };

	std::cout << MPField_MFEM << std::endl;

	EXPECT_TRUE(MPField_MFEM.isApprox(MPField_Hesthaven, tol_));

}


