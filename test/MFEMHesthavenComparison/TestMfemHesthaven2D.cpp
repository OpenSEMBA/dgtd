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

	SparseMatrix operatorToSparseMatrix(const Operator* op)
	{

		int width = op->Width();
		int height = op->Height();
		SparseMatrix res(height, height);
		Vector x(width), y(height);

		x = 0.0;

		for (int i = 0; i < width; i++)
		{
			x(i) = 1.0;
			op->Mult(x, y);
			for (int j = 0; j < height; j++)
			{
				if (y(j) != 0.0)
				{
					res.Add(i, j, y[j]);
				}
			}
			x(i) = 0.0;
		}

		res.Finalize();
		return res;
	}

	std::unique_ptr<SparseMatrix> rotateMatrixLexico(BilinearForm& matrix)
	{
		const Operator* rotatorOperator = matrix.FESpace()->GetElementRestriction(ElementDofOrdering::LEXICOGRAPHIC);
		const SparseMatrix rotatorMatrix = operatorToSparseMatrix(rotatorOperator);
		const SparseMatrix matrixSparse = matrix.SpMat();
		SparseMatrix* res;
		{
			auto aux = Mult(matrixSparse, rotatorMatrix);
			res = TransposeMult(rotatorMatrix, *aux);
		}
		return std::unique_ptr<SparseMatrix>(res);
	}

	Eigen::Matrix<double, 6, 6> rotatorO2{
	{0,0,0,0,0,1},
	{0,0,0,1,0,0},
	{1,0,0,0,0,0},
	{0,0,0,0,1,0},
	{0,1,0,0,0,0},
	{0,0,1,0,0,0}
	};
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

	EXPECT_TRUE(MFEMMass.isApprox(hesthavenMass * scaleFactor,tol_));
}

TEST_F(MFEMHesthaven2D, DOperators2D)
{
	setFES(2);

	Eigen::MatrixXd DrOperatorMFEM{
		buildInverseMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(X,*fes_)
	};

	Eigen::MatrixXd DrOperatorHesthaven{
		{ -1.5,  2.0, -0.5,  0.0, 0.0, 0.0},
		{ -0.5,  0.0,  0.5,  0.0, 0.0, 0.0},
		{  0.5, -2.0,  1.5,  0.0, 0.0, 0.0},
		{ -0.5,  1.0, -0.5, -1.0, 1.0, 0.0},
		{  0.5, -1.0,  0.5, -1.0, 1.0, 0.0},
		{  0.5,  0.0, -0.5, -2.0, 2.0, 0.0}
	};


	Eigen::Matrix<double, 6, 6> rotatedDrHesthaven = rotatorO2.transpose() * DrOperatorHesthaven * rotatorO2;

	std::cout << rotatedDrHesthaven << std::endl;

	Eigen::MatrixXd DsOperatorMFEM{
		buildInverseMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(Y,*fes_)
	};

	Eigen::MatrixXd DsOperatorHesthaven{
		{ -1.5,  0.0, 0.0,  2.0, 0.0, -0.5},
		{ -0.5, -1.0, 0.0,  1.0, 1.0, -0.5},
		{  0.5, -2.0, 0.0,  0.0, 2.0, -0.5},
		{ -0.5,  0.0, 0.0, -0.5, 0.0,  0.5},
		{  0.5, -1.0, 0.0, -1.0, 1.0,  0.5},
		{  0.5,  0.0, 0.0, -2.0, 0.0,  1.5}
	};

	Eigen::Matrix<double, 6, 6> rotatedDsHesthaven = rotatorO2.transpose() * DsOperatorHesthaven * rotatorO2;

	std::cout << rotatedDsHesthaven << std::endl;

	const double scaleFactor = 0.25;

	std::cout << DrOperatorMFEM << std::endl;
	std::cout << DsOperatorMFEM << std::endl;

	EXPECT_TRUE(DrOperatorMFEM.isApprox(scaleFactor * rotatedDrHesthaven, tol_));
	EXPECT_TRUE(DsOperatorMFEM.isApprox(scaleFactor * DsOperatorHesthaven, tol_));

}

TEST_F(MFEMHesthaven2D, manualMeshComparison)
{
	Mesh meshManual = Mesh::LoadFromFile("./TestData/twotriang.mesh", 1, 1);
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	Mesh meshAuto = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE, true);
	std::unique_ptr<FiniteElementCollection> fecAuto = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesAuto = std::make_unique<FiniteElementSpace>(&meshAuto, fecAuto.get());

	ASSERT_TRUE(buildMassMatrixEigen(*fesManual).isApprox(buildMassMatrixEigen(*fesAuto),tol_));
	ASSERT_TRUE(buildNormalStiffnessMatrixEigen(X, *fesManual).isApprox(buildNormalStiffnessMatrixEigen(X, *fesAuto), tol_));
	ASSERT_TRUE(buildNormalStiffnessMatrixEigen(Y, *fesManual).isApprox(buildNormalStiffnessMatrixEigen(Y, *fesAuto), tol_));
	
}

TEST_F(MFEMHesthaven2D, nodalPosition)
{
	Mesh meshAuto = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE, true);
	std::unique_ptr<FiniteElementCollection> fecAuto = std::make_unique<DG_FECollection>(2, 2, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesAuto = std::make_unique<FiniteElementSpace>(&meshAuto, fecAuto.get(), 2, Ordering::byNODES);

	GridFunction mfemNodes(fesAuto.get());
	meshAuto.GetNodes(mfemNodes);
	
	Eigen::Matrix<double, 6, 6> rotator{
		{0,0,0,0,0,1},
		{0,0,0,1,0,0},
		{1,0,0,0,0,0},
		{0,0,0,0,1,0},
		{0,1,0,0,0,0},
		{0,0,1,0,0,0}
	};
	Eigen::Matrix<double, 6, 6> identity;
	identity.setIdentity();

	Eigen::Matrix<double, 24, 24> fullRotator;
	fullRotator.setZero();
	fullRotator.block( 0,  0, 6, 6) = rotator;
	fullRotator.block( 6,  6, 6, 6) = identity;
	fullRotator.block(12, 12, 6, 6) = rotator;
	fullRotator.block(18, 18, 6, 6) = identity;

	Eigen::Vector<double, 24> mfemNodesPreRot;
	for (int i = 0; i < mfemNodes.Size(); i++) {
		mfemNodesPreRot(i, 0) = mfemNodes.Elem(i);
	}

	Eigen::Vector<double, 24> rotatedMfemNodesVector = fullRotator * mfemNodesPreRot;

	Eigen::Matrix<double, 12, 2> rotatedMfemNodes;
	for (int i = 0; i < mfemNodes.Size()/2; i++) {
		rotatedMfemNodes(i, 0) = rotatedMfemNodesVector(i);
		rotatedMfemNodes(i, 1) = rotatedMfemNodesVector(i + rotatedMfemNodesVector.rows()/2);
	}

	Eigen::MatrixXd hesthavenNodes{
		{ 0.0, 1.0},
		{ 0.0, 0.5},
		{ 0.0, 0.0},
		{ 0.5, 1.0},
		{ 0.5, 0.5},
		{ 1.0, 1.0},
		{ 1.0, 1.0},
		{ 0.5, 0.5},
		{ 0.0, 0.0},
		{ 1.0, 0.5},
		{ 0.5, 0.0},
		{ 1.0, 0.0}
	};

	EXPECT_TRUE(hesthavenNodes.isApprox(rotatedMfemNodes));
}
