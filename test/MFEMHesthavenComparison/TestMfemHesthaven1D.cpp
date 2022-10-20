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

TEST_F(MFEMHesthaven1D, MSOperator)
{
	setFES(2, 4);
	auto MS_MFEM4E = toEigen(*buildByMult(
		*buildInverseMassMatrix(E, Model(mesh_, AttributeToMaterial{}, { {1, BdrCond::SMA}, {2, BdrCond::SMA} }), *fes_),
		*buildDerivativeOperator(X, *fes_), *fes_)
		.get()->SpMat().ToDenseMatrix());
	auto MS_Hesthaven4E = buildMatrixForMSTest4E();

	setFES(2, 3);
	auto MS_MFEM3E = toEigen(*buildByMult(
		*buildInverseMassMatrix(E, Model(mesh_, AttributeToMaterial{}, { {1, BdrCond::SMA}, {2, BdrCond::SMA} }), *fes_),
		*buildDerivativeOperator(X, *fes_), *fes_)
		.get()->SpMat().ToDenseMatrix());
	Eigen::MatrixXd MS_Hesthaven3E{
		{  9.0, -12.0,  3.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0},
		{  3.0,   0.0, -3.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0},
		{ -3.0,  12.0, -9.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0},
		{  0.0,   0.0,  0.0,  9.0, -12.0,  3.0,  0.0,   0.0,  0.0},
		{  0.0,   0.0,  0.0,  3.0,   0.0, -3.0,  0.0,   0.0,  0.0},
		{  0.0,   0.0,  0.0, -3.0,  12.0, -9.0,  0.0,   0.0,  0.0},
		{  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  9.0, -12.0,  3.0},
		{  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  3.0,   0.0, -3.0},
		{  0.0,   0.0,  0.0,  0.0,   0.0,  0.0, -3.0,  12.0, -9.0},
	};

	double multMS_coefficient = -1.0;

	EXPECT_TRUE(MS_Hesthaven4E.isApprox(multMS_coefficient * MS_MFEM4E, tol_));
	EXPECT_TRUE(MS_Hesthaven3E.isApprox(multMS_coefficient * MS_MFEM3E, tol_));
}

TEST_F(MFEMHesthaven1D, MFOperatorPEC_E) //TODO
{
	setFES(2, 4);
	Eigen::MatrixXd MF_MFEM4E{
		buildInverseMassMatrixEigen(*fes_) * buildNormalPECFluxOperator1D(E, *fes_, std::vector<int>{0})
	};

	Eigen::MatrixXd MF_Hesthaven4E{
		{-36.0, 0.0,  6.0,  -6.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  6.0, 0.0, -3.0,   3.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{-12.0, 0.0, 18.0, -18.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, 18.0, -18.0, 0.0,  6.0,  -6.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -3.0,   3.0, 0.0, -3.0,   3.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  6.0,  -6.0, 0.0, 18.0, -18.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, 18.0, -18.0, 0.0,   6.0,  -6.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, -3.0,   3.0, 0.0,  -3.0,   3.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  6.0,  -6.0, 0.0,  18.0, -18.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  18.0, -18.0, 0.0,  12.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,  -6.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,   6.0,  -6.0, 0.0,  36.0}
	};

	std::cout << MF_MFEM4E << std::endl;

	setFES(2, 3);
	Eigen::MatrixXd MF_MFEM3E{
		buildInverseMassMatrixEigen(*fes_) * buildNormalPECFluxOperator1D(E, *fes_, std::vector<int>{0})
	};

	Eigen::MatrixXd MF_Hesthaven3E{
		{-27.0, 0.0,   4.5,  -4.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  4.5, 0.0, -2.25,  2.25, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ -9.0, 0.0,  13.5, -13.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  13.5, -13.5, 0.0,   4.5,  -4.5, 0.0,   0.0},
		{  0.0, 0.0, -2.25,  2.25, 0.0, -2.25,  2.25, 0.0,   0.0},
		{  0.0, 0.0,   4.5,  -4.5, 0.0,  13.5, -13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  13.5, -13.5, 0.0,   9.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0, -2.25,  2.25, 0.0,  -4.5},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   4.5,  -4.5, 0.0,  27.0}
	};

	EXPECT_TRUE(MF_MFEM4E.isApprox(MF_Hesthaven4E, tol_));
	EXPECT_TRUE(MF_MFEM3E.isApprox(MF_Hesthaven3E, tol_));
}

TEST_F(MFEMHesthaven1D, MFOperatorPEC_H) //TODO
{
	setFES(2, 4);
	Eigen::MatrixXd MF_MFEM4E{
		buildInverseMassMatrixEigen(*fes_) * buildNormalPECFluxOperator1D(H, *fes_, std::vector<int>{0})
	};

	Eigen::MatrixXd MF_Hesthaven4E{
		{  0.0, 0.0,  6.0,  -6.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -3.0,   3.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, 18.0, -18.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, 18.0, -18.0, 0.0,  6.0,  -6.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -3.0,   3.0, 0.0, -3.0,   3.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  6.0,  -6.0, 0.0, 18.0, -18.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, 18.0, -18.0, 0.0,   6.0,  -6.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, -3.0,   3.0, 0.0,  -3.0,   3.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  6.0,  -6.0, 0.0,  18.0, -18.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  18.0, -18.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,   6.0,  -6.0, 0.0,   0.0}
	};

	setFES(2, 3);
	Eigen::MatrixXd MF_MFEM3E{
		buildInverseMassMatrixEigen(*fes_) * buildNormalPECFluxOperator1D(H, *fes_, std::vector<int>{0})
	};

	Eigen::MatrixXd MF_Hesthaven3E{
		{  0.0, 0.0,   4.5,  -4.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -2.25,  2.25, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  13.5, -13.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  13.5, -13.5, 0.0,   4.5,  -4.5, 0.0,   0.0},
		{  0.0, 0.0, -2.25,  2.25, 0.0, -2.25,  2.25, 0.0,   0.0},
		{  0.0, 0.0,   4.5,  -4.5, 0.0,  13.5, -13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  13.5, -13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0, -2.25,  2.25, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   4.5,  -4.5, 0.0,   0.0}
	};

	EXPECT_TRUE(MF_MFEM4E.isApprox(MF_Hesthaven4E, tol_));
	EXPECT_TRUE(MF_MFEM3E.isApprox(MF_Hesthaven3E, tol_));
}

TEST_F(MFEMHesthaven1D, MFOperatorSMA)
{
	setFES(2, 4);
	Eigen::MatrixXd MF_MFEM4E{
		buildInverseMassMatrixEigen(*fes_) * buildNormalSMAFluxOperator1D(*fes_, std::vector<int>{0})
	};

	Eigen::MatrixXd MF_Hesthaven4E{
		{-18.0, 0.0,  6.0,  -6.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  3.0, 0.0, -3.0,   3.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ -6.0, 0.0, 18.0, -18.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, 18.0, -18.0, 0.0,  6.0,  -6.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -3.0,   3.0, 0.0, -3.0,   3.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  6.0,  -6.0, 0.0, 18.0, -18.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, 18.0, -18.0, 0.0,   6.0,  -6.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, -3.0,   3.0, 0.0,  -3.0,   3.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  6.0,  -6.0, 0.0,  18.0, -18.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  18.0, -18.0, 0.0,   6.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,  -3.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,   6.0,  -6.0, 0.0,  18.0}
	};

	setFES(2, 3);
	Eigen::MatrixXd MF_MFEM3E{
		buildInverseMassMatrixEigen(*fes_) * buildNormalSMAFluxOperator1D(*fes_, std::vector<int>{0})
	};

	Eigen::MatrixXd MF_Hesthaven3E{
		{-13.5, 0.0,   4.5,  -4.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ 2.25, 0.0, -2.25,  2.25, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ -4.5, 0.0,  13.5, -13.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  13.5, -13.5, 0.0,   4.5,  -4.5, 0.0,   0.0},
		{  0.0, 0.0, -2.25,  2.25, 0.0, -2.25,  2.25, 0.0,   0.0},
		{  0.0, 0.0,   4.5,  -4.5, 0.0,  13.5, -13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  13.5, -13.5, 0.0,   4.5},
		{  0.0, 0.0,   0.0,   0.0, 0.0, -2.25,  2.25, 0.0, -2.25},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   4.5,  -4.5, 0.0,  13.5}
	};

	EXPECT_TRUE(MF_MFEM4E.isApprox(MF_Hesthaven4E, tol_));
	EXPECT_TRUE(MF_MFEM3E.isApprox(MF_Hesthaven3E, tol_));
} 


TEST_F(MFEMHesthaven1D, MPOperatorPEC_E) //TODO
{
	setFES(2, 4);
	Eigen::MatrixXd MP_MFEM4E{
		buildInverseMassMatrixEigen(*fes_) * buildPECPenaltyOperator1D(E, *fes_)
	};

	Eigen::MatrixXd MP_Hesthaven4E{
		{-36.0, 0.0, -6.0,   6.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  6.0, 0.0,  3.0,  -3.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{-12.0, 0.0,-18.0,  18.0, 0.0,  0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, 18.0, -18.0, 0.0, -6.0,   6.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -3.0,   3.0, 0.0,  3.0,  -3.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  6.0,  -6.0, 0.0,-18.0,  18.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, 18.0, -18.0, 0.0,  -6.0,   6.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0, -3.0,   3.0, 0.0,   3.0,  -3.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  6.0,  -6.0, 0.0, -18.0,  18.0, 0.0,   0.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  18.0, -18.0, 0.0, -12.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,   6.0},
		{  0.0, 0.0,  0.0,   0.0, 0.0,  0.0,   0.0, 0.0,   6.0,  -6.0, 0.0, -36.0}
	};

	std::cout << MP_MFEM4E << std::endl;

	setFES(2, 3);
	Eigen::MatrixXd MP_MFEM3E{
		buildInverseMassMatrixEigen(*fes_)* buildPECPenaltyOperator1D(E, *fes_)
	};

	Eigen::MatrixXd MP_Hesthaven3E{
		{-27.0, 0.0,  -4.5,   4.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  4.5, 0.0,  2.25, -2.25, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ -9.0, 0.0, -13.5,  13.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  13.5, -13.5, 0.0,  -4.5,   4.5, 0.0,   0.0},
		{  0.0, 0.0, -2.25,  2.25, 0.0,  2.25, -2.25, 0.0,   0.0},
		{  0.0, 0.0,   4.5,  -4.5, 0.0, -13.5,  13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  13.5, -13.5, 0.0,  -9.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0, -2.25,  2.25, 0.0,   4.5},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   4.5,  -4.5, 0.0, -27.0}
	};

	double multMP_coefficient = -1.0;

	EXPECT_TRUE(MP_MFEM4E.isApprox(multMP_coefficient * MP_Hesthaven4E, tol_));
	EXPECT_TRUE(MP_MFEM3E.isApprox(multMP_coefficient * MP_Hesthaven3E, tol_));
}


TEST_F(MFEMHesthaven1D, MPOperatorPEC_H) //TODO
{
	setFES(2, 4);
	Eigen::MatrixXd MP_MFEM4E{
		buildInverseMassMatrixEigen(*fes_) * buildPECPenaltyOperator1D(H, *fes_)
	};

	Eigen::MatrixXd MP_Hesthaven4E{
		{  0.0, 0.0,  -6.0,   6.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,   3.0,  -3.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -18.0,  18.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  18.0, -18.0, 0.0,  -6.0,   6.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  -3.0,   3.0, 0.0,   3.0,  -3.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,   6.0,  -6.0, 0.0, -18.0,  18.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  18.0, -18.0, 0.0,  -6.0,   6.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,   3.0,  -3.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   6.0,  -6.0, 0.0, -18.0,  18.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,  18.0, -18.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   6.0,  -6.0, 0.0,   0.0}
	};

	setFES(2, 3);
	Eigen::MatrixXd MP_MFEM3E{
		buildInverseMassMatrixEigen(*fes_) * buildPECPenaltyOperator1D(H, *fes_)
	};

	Eigen::MatrixXd MP_Hesthaven3E{
		{  0.0, 0.0,  -4.5,   4.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  2.25, -2.25, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0, -13.5,  13.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  13.5, -13.5, 0.0,  -4.5,   4.5, 0.0,   0.0},
		{  0.0, 0.0, -2.25,  2.25, 0.0,  2.25, -2.25, 0.0,   0.0},
		{  0.0, 0.0,   4.5,  -4.5, 0.0, -13.5,  13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  13.5, -13.5, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0, -2.25,  2.25, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   4.5,  -4.5, 0.0,   0.0}
	};

	double multMP_coefficient = -1.0;

	EXPECT_TRUE(MP_MFEM4E.isApprox(multMP_coefficient * MP_Hesthaven4E, tol_));
	EXPECT_TRUE(MP_MFEM3E.isApprox(multMP_coefficient * MP_Hesthaven3E, tol_));
}
TEST_F(MFEMHesthaven1D, MPOperatorSMA)
{
	setFES(2, 4);
	Eigen::MatrixXd MP_MFEM4E{
		buildInverseMassMatrixEigen(*fes_) * buildSMAPenaltyOperator1D(*fes_)
	};
	Eigen::MatrixXd MP_Hesthaven4E{
		{-18.0, 0.0,  -6.0,   6.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  3.0, 0.0,   3.0,  -3.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ -6.0, 0.0, -18.0,  18.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  18.0, -18.0, 0.0,  -6.0,   6.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,  -3.0,   3.0, 0.0,   3.0,  -3.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,   6.0,  -6.0, 0.0, -18.0,  18.0, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  18.0, -18.0, 0.0,  -6.0,   6.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,   3.0,  -3.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   6.0,  -6.0, 0.0, -18.0,  18.0, 0.0,   0.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,  18.0, -18.0, 0.0,  -6.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,  -3.0,   3.0, 0.0,   3.0},
		{  0.0, 0.0,   0.0,   0.0, 0.0,   0.0,   0.0, 0.0,   6.0,  -6.0, 0.0, -18.0}
	};

	setFES(2, 3);
	Eigen::MatrixXd MP_MFEM3E{
		buildInverseMassMatrixEigen(*fes_) * buildSMAPenaltyOperator1D(*fes_)
	};

	Eigen::MatrixXd MP_Hesthaven3E{
		{-13.5, 0.0,   -4.5,   4.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ 2.25, 0.0,   2.25, -2.25, 0.0,   0.0,   0.0, 0.0,   0.0},
		{ -4.5, 0.0,  -13.5,  13.5, 0.0,   0.0,   0.0, 0.0,   0.0},
		{  0.0, 0.0,   13.5, -13.5, 0.0,  -4.5,   4.5, 0.0,   0.0},
		{  0.0, 0.0,  -2.25,  2.25, 0.0,  2.25, -2.25, 0.0,   0.0},
		{  0.0, 0.0,    4.5,  -4.5, 0.0, -13.5,  13.5, 0.0,   0.0},
		{  0.0, 0.0,    0.0,   0.0, 0.0,  13.5, -13.5, 0.0,  -4.5},
		{  0.0, 0.0,    0.0,   0.0, 0.0, -2.25,  2.25, 0.0,  2.25},
		{  0.0, 0.0,    0.0,   0.0, 0.0,   4.5,  -4.5, 0.0, -13.5}
	};

	double multMP_coefficient = -1.0;

	EXPECT_TRUE(MP_Hesthaven4E.isApprox(multMP_coefficient * MP_MFEM4E, tol_));
	EXPECT_TRUE(MP_Hesthaven3E.isApprox(multMP_coefficient * MP_MFEM3E, tol_));
}
