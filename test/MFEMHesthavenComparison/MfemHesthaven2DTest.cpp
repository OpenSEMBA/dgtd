#include <gtest/gtest.h>

#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/Types.h"
#include "MfemHesthavenFunctionsTest.h"
#include "GlobalFunctions.h"
#include "AnalyticalFunctions2D.h"
#include "SourceFixtures.h"
#include "maxwell/Solver.h"

using namespace mfem;
using namespace maxwell;
using namespace maxwell::fixtures::sources;

class MFEMHesthaven2D : public ::testing::Test {
protected:

	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE);
		fec_ = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get(), 1, 0);
	}

	void setFES(
		const int order, 
		const int nx = 1, const int ny = 1, 
		const double sx = 1.0, const double sy = 1.0, 
		const int vdim = 1)
	{
		mesh_ = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, false, sx, sy);
		fec_ = std::make_unique<DG_FECollection>(order, 2, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get(), vdim, 0);
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

	Probes buildExportProbes()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}
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
	std::cout << hesthavenMass * scaleFactor << std::endl;

	EXPECT_TRUE(MFEMMass.isApprox(hesthavenMass * scaleFactor,1e-4));
}

TEST_F(MFEMHesthaven2D, DrOperator2D)
{
	setFES(2);

	Eigen::MatrixXd DrOperatorHesthaven{
	{ -1.5,  2.0, -0.5,  0.0, 0.0, 0.0},
	{ -0.5,  0.0,  0.5,  0.0, 0.0, 0.0},
	{  0.5, -2.0,  1.5,  0.0, 0.0, 0.0},
	{ -0.5,  1.0, -0.5, -1.0, 1.0, 0.0},
	{  0.5, -1.0,  0.5, -1.0, 1.0, 0.0},
	{  0.5,  0.0, -0.5, -2.0, 2.0, 0.0}
	};

	Eigen::Matrix<double, 6, 6> rotatedDrHesthaven = rotatorO2.transpose() * DrOperatorHesthaven * rotatorO2;
	Eigen::Matrix<double, 12, 12> globalDrHesthaven;
	globalDrHesthaven.setZero();
	globalDrHesthaven.block(0, 0, 6, 6) = -1.0 * rotatedDrHesthaven;
	globalDrHesthaven.block(6, 6, 6, 6) = rotatedDrHesthaven;

	Eigen::MatrixXd DrOperatorMFEM{
		buildInverseMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(Y,*fes_)
	};

	std::cout << globalDrHesthaven << std::endl;
	std::cout << DrOperatorMFEM << std::endl;

	EXPECT_TRUE(DrOperatorMFEM.isApprox(globalDrHesthaven, tol_));

}
TEST_F(MFEMHesthaven2D, DsOperator2D)
{
	setFES(2);

	Eigen::MatrixXd DsOperatorHesthaven{
	{ -1.5,  0.0, 0.0,  2.0, 0.0, -0.5},
	{ -0.5, -1.0, 0.0,  1.0, 1.0, -0.5},
	{  0.5, -2.0, 0.0,  0.0, 2.0, -0.5},
	{ -0.5,  0.0, 0.0,  0.0, 0.0,  0.5},
	{  0.5, -1.0, 0.0, -1.0, 1.0,  0.5},
	{  0.5,  0.0, 0.0, -2.0, 0.0,  1.5}
	};

	Eigen::Matrix<double, 6, 6> rotatedDsHesthaven = rotatorO2.transpose() * DsOperatorHesthaven * rotatorO2;
	Eigen::Matrix<double, 12, 12> globalDsHesthaven;
	globalDsHesthaven.setZero();
	globalDsHesthaven.block(0, 0, 6, 6) = rotatedDsHesthaven;
	globalDsHesthaven.block(6, 6, 6, 6) = -1.0 * rotatedDsHesthaven;

	Eigen::MatrixXd DsOperatorMFEM{
		buildInverseMassMatrixEigen(*fes_) * buildNormalStiffnessMatrixEigen(X,*fes_)
	};

	std::cout << globalDsHesthaven << std::endl;
	std::cout << DsOperatorMFEM << std::endl;

	EXPECT_TRUE(DsOperatorMFEM.isApprox(globalDsHesthaven, tol_));

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

	std::cout << mfemNodes << std::endl;

	Eigen::Matrix<double, 6, 6> identity;
	identity.setIdentity();

	Eigen::Matrix<double, 24, 24> fullRotator;
	fullRotator.setZero();
	fullRotator.block( 0,  0, 6, 6) = rotatorO2;
	fullRotator.block( 6,  6, 6, 6) = identity;
	fullRotator.block(12, 12, 6, 6) = rotatorO2;
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
TEST_F(MFEMHesthaven2D, oneFace)
{
	Mesh meshManual = Mesh::LoadFromFile("./TestData/onetriang.mesh", true, 1);
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
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
		Mesh meshOne = Mesh::MakeCartesian1D(1);
		std::unique_ptr<FiniteElementCollection> fecOne = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		std::unique_ptr<FiniteElementSpace> fesOne = std::make_unique<FiniteElementSpace>(&meshOne, fecOne.get());

	}

}

TEST_F(MFEMHesthaven2D, DISABLED_MFEMHesthavenSameMesh)
{
	Mesh meshManual = Mesh::LoadFromFile("./TestData/2dplane.msh", true, 1);
	std::unique_ptr<FiniteElementCollection> fecManual = std::make_unique<DG_FECollection>(4, 2, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fesManual = std::make_unique<FiniteElementSpace>(&meshManual, fecManual.get());

	Model model = Model(meshManual, AttributeToMaterial{}, AttributeToBoundary{});

	maxwell::Solver solver{
		model,
		buildExportProbes(),
		buildGaussianInitialField(),
		SolverOptions{}
			.setTimeStep(0.012587)
			.setFinalTime(1.0)
			.setOrder(10)
	};

}