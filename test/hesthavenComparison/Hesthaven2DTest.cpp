#include <gtest/gtest.h>

#include "TestUtils.h"
#include "HesthavenFunctions.h"
#include "math/EigenMfemTools.h"
#include "evolution/EvolutionMethods.h"

using namespace mfem;
using namespace maxwell;

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

	Eigen::Matrix<double, 6, 6> rotatorO2 {
		{0,0,0,0,0,1},
		{0,0,0,1,0,0},
		{1,0,0,0,0,0},
		{0,0,0,0,1,0},
		{0,1,0,0,0,0},
		{0,0,1,0,0,0}
	};

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}

	std::pair<FiniteElementSpace, Model> buildRequirementsForComparison()
	{
		auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
		auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
		auto fes{ FiniteElementSpace(&mesh,&fec) };

		std::pair<FiniteElementSpace, Model> res(
			fes, 
			Model(
				mesh,
				GeomTagToMaterial(),
				GeomTagToBoundary(),
				GeomTagToInteriorConditions())
		);
		return res;
	}
};

TEST_F(MFEMHesthaven2D, massMatrix)
{
	// Hesthaven's mass matrix is calculated with
	// $$ Mass = [V V']^{-1} J $$
	// where $V$ is the Vandermonde matrix and $J$ is the jacobian.
	// His reference element has area $A_r = 2$.
	// In this FiniteElementSpace there are two triangles from
	// [0, 1] x [0, 1]. Therefore they have A = 1/2.
	// The jacobian is 
	// $$ J = A / A_r = 0.25 $$

	Eigen::MatrixXd vanderProdInverse{
		{0.33333, 0.16667, 0.16667,     0.0,     0.0,     0.0},
		{0.16667, 0.33333, 0.16667,     0.0,     0.0,     0.0},
		{0.16667, 0.16667, 0.33333,     0.0,     0.0,     0.0},
		{    0.0,     0.0,     0.0, 0.33333, 0.16667, 0.16667},
		{    0.0,     0.0,     0.0, 0.16667, 0.33333, 0.16667},
		{    0.0,     0.0,     0.0, 0.16667, 0.16667, 0.33333}
	};
	const double jacobian = 0.25;

	auto hesthavenMass{ vanderProdInverse * jacobian };
	auto MFEMMass = buildMassMatrixEigen(*fes_);

	EXPECT_TRUE(MFEMMass.isApprox(hesthavenMass, 1e-4));
}

TEST_F(MFEMHesthaven2D, DrOperator)
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

	EXPECT_TRUE(DrOperatorMFEM.isApprox(globalDrHesthaven, tol_));

}

TEST_F(MFEMHesthaven2D, 2D_Operator_ZeroNormal_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh,GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd ZeroNormalOperator{
		{ 10., -1.12132034, -1.12132034,  2.12132034,  2.12132034,  0.},
		{ -2.,  8.53553391, -2.29289322, -0.70710678, -3.53553391,  0.},
		{ -2., -2.29289322,  8.53553391, -3.53553391, -0.70710678,  0.},
		{  0., -0.70710678, -3.53553391,  8.53553391, -2.29289322, -2.},
		{  0., -3.53553391, -0.70710678, -2.29289322,  8.53553391, -2.},
		{  0.,  2.12132034,  2.12132034, -1.12132034, -1.12132034, 10.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMP = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(E, model, fes),
			*buildZeroNormalOperator(E, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	for (int i = 0; i < EigenMP.rows(); i++) {
		for (int j = 0; j < EigenMP.cols(); j++) {
			ASSERT_TRUE(abs(EigenMP(i, j) - ZeroNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_OneNormal_nxEZ_HX_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd OneNormalOperator{
		{ 0., -1.5, -1.5,  1.5,  1.5, 0.},
		{ 0.,  2.5,  0.5, -0.5, -2.5, 0.},
		{ 0.,  0.5,  2.5, -2.5, -0.5, 0.},
		{ 0.,  0.5,  2.5, -2.5, -0.5, 0.},
		{ 0.,  2.5,  0.5, -0.5, -2.5, 0.},
		{ 0., -1.5, -1.5,  1.5,  1.5, 0.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(E, model, fes),
			*buildOneNormalOperator(H, { X }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFN << std::endl;

	for (int i = 0; i < EigenMFN.rows(); i++) {
		for (int j = 0; j < EigenMFN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFN(i, j) - OneNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_OneNormal_nyEZ_HY_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd OneNormalOperator{
		{ 0.,  1.5,  1.5, -1.5, -1.5, 0.},
		{ 0., -2.5, -0.5,  0.5,  2.5, 0.},
		{ 0., -0.5, -2.5,  2.5,  0.5, 0.},
		{ 0., -0.5, -2.5,  2.5,  0.5, 0.},
		{ 0., -2.5, -0.5,  0.5,  2.5, 0.},
		{ 0.,  1.5,  1.5, -1.5, -1.5, 0.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(E, model, fes),
			*buildOneNormalOperator(H, { Y }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFN << std::endl;

	for (int i = 0; i < EigenMFN.rows(); i++) {
		for (int j = 0; j < EigenMFN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFN(i, j) - OneNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_OneNormal_nyHX_EZ_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd OneNormalOperator{
		{ 5.,  1.5,  2.5, -1.5, -1.5,  0.},
		{-3., -2.5, -3.5,  0.5,  2.5,  0.},
		{ 1., -0.5,  2.5,  2.5,  0.5,  0.},
		{ 0., -0.5, -2.5,  2.5,  3.5,  3.},
		{ 0., -2.5, -0.5,  0.5, -2.5, -1.},
		{ 0.,  1.5,  1.5, -1.5, -2.5, -5.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(H, model, fes),
			*buildOneNormalOperator(E, { Y }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFN << std::endl;

	for (int i = 0; i < EigenMFN.rows(); i++) {
		for (int j = 0; j < EigenMFN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFN(i, j) - OneNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_OneNormal_nxHY_EZ_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd OneNormalOperator{
		{-5., -2.5, -1.5,  1.5,  1.5,  0.},
		{-1., -2.5,  0.5, -0.5, -2.5,  0.},
		{ 3.,  3.5,  2.5, -2.5, -0.5,  0.},
		{ 0.,  0.5,  2.5,  2.5, -0.5,  1.},
		{ 0.,  2.5,  0.5, -3.5, -2.5, -3.},
		{ 0., -1.5, -1.5,  2.5,  1.5,  5.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(H, model, fes),
			*buildOneNormalOperator(E, { X }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFN << std::endl;

	for (int i = 0; i < EigenMFN.rows(); i++) {
		for (int j = 0; j < EigenMFN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFN(i, j) - OneNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_TwoNormal_nxHXnx_HX_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd TwoNormalOperator{
		{ 0., -1.06066017, -1.06066017,  1.06066017,  1.06066017,  0.},
		{ 0.,  1.76776695,  0.35355339, -0.35355339, -1.76776695,  0.},
		{ 0.,  0.35355339,  1.76776695, -1.76776695, -0.35355339,  0.},
		{ 0., -0.35355339, -1.76776695,  1.76776695,  0.35355339,  0.},
		{ 0., -1.76776695, -0.35355339,  0.35355339,  1.76776695,  0.},
		{ 0.,  1.06066017,  1.06066017, -1.06066017, -1.06066017,  0.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFNN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(H, model, fes),
			*buildTwoNormalOperator(H, { X, X }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFNN << std::endl;

	for (int i = 0; i < EigenMFNN.rows(); i++) {
		for (int j = 0; j < EigenMFNN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFNN(i, j) - TwoNormalOperator(i, j)) < 1e-3);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_TwoNormal_nxHXny_HY_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd TwoNormalOperator{
		{ 0.,  1.06066017,  1.06066017, -1.06066017, -1.06066017, 0.},
		{ 0., -1.76776695, -0.35355339,  0.35355339,  1.76776695, 0.},
		{ 0., -0.35355339, -1.76776695,  1.76776695,  0.35355339, 0.},
		{ 0.,  0.35355339,  1.76776695, -1.76776695, -0.35355339, 0.},
		{ 0.,  1.76776695,  0.35355339, -0.35355339, -1.76776695, 0.},
		{ 0., -1.06066017, -1.06066017,  1.06066017,  1.06066017, 0.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFNN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(H, model, fes),
			*buildTwoNormalOperator(H, { X, Y }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFNN << std::endl;

	for (int i = 0; i < EigenMFNN.rows(); i++) {
		for (int j = 0; j < EigenMFNN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFNN(i, j) - TwoNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_TwoNormal_nyHYnx_HY_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd TwoNormalOperator{
		{ 0.,  1.06066017,  1.06066017, -1.06066017, -1.06066017, 0.},
		{ 0., -1.76776695, -0.35355339,  0.35355339,  1.76776695, 0.},
		{ 0., -0.35355339, -1.76776695,  1.76776695,  0.35355339, 0.},
		{ 0.,  0.35355339,  1.76776695, -1.76776695, -0.35355339, 0.},
		{ 0.,  1.76776695,  0.35355339, -0.35355339, -1.76776695, 0.},
		{ 0., -1.06066017, -1.06066017,  1.06066017,  1.06066017, 0.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFNN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(H, model, fes),
			*buildTwoNormalOperator(H, { Y, X }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFNN << std::endl;

	for (int i = 0; i < EigenMFNN.rows(); i++) {
		for (int j = 0; j < EigenMFNN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFNN(i, j) - TwoNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, 2D_Operator_TwoNormal_nyHYny_HY_PEC)
{

	auto mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "Maxwell2D_K2.mesh").c_str(),1,0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	GeomTagToBoundary pecBdr{ {2,BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterial(), pecBdr, GeomTagToInteriorConditions());

	Eigen::MatrixXd TwoNormalOperator{
		{ 0., -1.06066017, -1.06066017,  1.06066017,  1.06066017,  0.},
		{ 0.,  1.76776695,  0.35355339, -0.35355339, -1.76776695,  0.},
		{ 0.,  0.35355339,  1.76776695, -1.76776695, -0.35355339,  0.},
		{ 0., -0.35355339, -1.76776695,  1.76776695,  0.35355339,  0.},
		{ 0., -1.76776695, -0.35355339,  0.35355339,  1.76776695,  0.},
		{ 0.,  1.06066017,  1.06066017, -1.06066017, -1.06066017,  0.}
	};

	EvolutionOptions opts = EvolutionOptions();
	opts.order = 1;
	auto EigenMFNN = toEigen(
		*buildByMult(
			*buildInverseMassMatrix(H, model, fes),
			*buildTwoNormalOperator(H, { Y, Y }, model, fes, opts),
			fes
		)->SpMat().ToDenseMatrix()
	);

	std::cout << EigenMFNN << std::endl;

	for (int i = 0; i < EigenMFNN.rows(); i++) {
		for (int j = 0; j < EigenMFNN.cols(); j++) {
			ASSERT_TRUE(abs(EigenMFNN(i, j) - TwoNormalOperator(i, j)) < 1e-8);
		}
	}

}

TEST_F(MFEMHesthaven2D, DsOperator)
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

	EXPECT_TRUE(DsOperatorMFEM.isApprox(globalDsHesthaven, tol_));

}

TEST_F(MFEMHesthaven2D, manualMeshComparison)
{
	Mesh meshManual = Mesh::LoadFromFile((mfemMeshes2DFolder() + "twotriang.mesh").c_str(), 1, 1);
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
		{ 0.0, 1.0 },
		{ 0.0, 0.5 },
		{ 0.0, 0.0 },
		{ 0.5, 1.0 },
		{ 0.5, 0.5 },
		{ 1.0, 1.0 },
		{ 1.0, 1.0 },
		{ 0.5, 0.5 },
		{ 0.0, 0.0 },
		{ 1.0, 0.5 },
		{ 0.5, 0.0 },
		{ 1.0, 0.0 }
	};

	EXPECT_TRUE(hesthavenNodes.isApprox(rotatedMfemNodes));
}