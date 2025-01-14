#include <gtest/gtest.h>

#include "TestUtils.h"
#include "HesthavenFunctions.h"
#include "components/Model.h"
#include "components/Types.h"
#include "evolution/EvolutionMethods.h"
#include "evolution/HesthavenEvolutionTools.h"
#include "math/EigenMfemTools.h"

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

	double tol_ = 1e-4;
	double hesthaven_triangle_scaling_factor = 0.25;

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
				GeomTagToMaterialInfo(),
				GeomTagToBoundaryInfo())
		);
		return res;
	}

	const double calculateElementJacobian(const FiniteElement* fe, ElementTransformation* e_trans) 
	{
		const int order = fe->GetOrder() + fe->GetOrder() + e_trans->OrderW();
		const IntegrationRule& ir = IntRules.Get(e_trans->GetGeometryType(),
			order);
		double res = 0.0;
		for (int i = 0; i < ir.GetNPoints(); i++)
		{
			const IntegrationPoint& ip = ir.IntPoint(i);
			e_trans->SetIntPoint(&ip);
			res += ip.weight * e_trans->Weight();
		}
		return res;
	}

	const int getFaceIntegrationOrder(const FaceElementTransformations* f_trans, const FiniteElementSpace& fes)
	{
		int res = f_trans->Elem1->OrderW() + 2 * fes.GetFE(f_trans->Elem1No)->GetOrder();
		if (fes.GetFE(f_trans->Elem1No)->Space() == FunctionSpace::Pk) {
			res++;
		}
		return res;
	}

	const double calculateSurfaceJacobian(const FiniteElement* fe, FaceElementTransformations* f_trans, const FiniteElementSpace& fes)
	{
		auto order = getFaceIntegrationOrder(f_trans, fes);
		const IntegrationRule& fir = IntRules.Get(f_trans->GetGeometryType(), order);
		double res = 0.0;
		for (int i = 0; i < fir.GetNPoints(); i++)
		{
			const IntegrationPoint& fip = fir.IntPoint(i);
			f_trans->SetIntPoint(&fip);
			res += fip.weight * f_trans->Weight();
		}
		return res;
	}

	DynamicMatrix assembleMassMatrix(FiniteElementSpace& fes)
	{
		BilinearForm bf(&fes);
		ConstantCoefficient one(1.0);
		bf.AddDomainIntegrator(new MassIntegrator(one));
		bf.Assemble();
		bf.Finalize();

		return toEigen(*bf.SpMat().ToDenseMatrix());
	}

	SubMesh assembleInteriorFaceSubMesh(Mesh& m_copy, const std::map<int, Attribute>& att_map)
	{
		Array<int> sm_tag;
		sm_tag.Append(hesthavenMeshingTag);
		auto faces = getFacesForElement(m_copy, 0);
		auto f_trans = getInteriorFaceTransformation(m_copy, faces);
		markElementsForSubMeshing(f_trans, m_copy);
		auto res = SubMesh::CreateFromDomain(m_copy, sm_tag);
		restoreOriginalAttributesAfterSubMeshing(f_trans, m_copy, att_map);
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

	DynamicMatrix vanderProdInverse{
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

	DynamicMatrix DrOperatorHesthaven{
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

	DynamicMatrix DrOperatorMFEM{
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
	Model model(mesh,GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix ZeroNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix OneNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix OneNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix OneNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix OneNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix TwoNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix TwoNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix TwoNormalOperator{
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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

	DynamicMatrix TwoNormalOperator{
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

	DynamicMatrix DsOperatorHesthaven{
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

	DynamicMatrix DsOperatorMFEM{
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

	DynamicMatrix hesthavenNodes{
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

TEST_F(MFEMHesthaven2D, triangleElementJacobianO1)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	ASSERT_NEAR(0.5, calculateElementJacobian(fes.GetFE(0), fes.GetElementTransformation(0)), tol_);
}

TEST_F(MFEMHesthaven2D, triangleElementJacobianO2)
{
	const int basis_order = 2;
	auto m{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	ASSERT_NEAR(0.5, calculateElementJacobian(fes.GetFE(0), fes.GetElementTransformation(0)), tol_);
}

TEST_F(MFEMHesthaven2D, triangleElementJacobianO3)
{
	const int basis_order = 3;
	auto m{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	ASSERT_NEAR(0.5, calculateElementJacobian(fes.GetFE(0), fes.GetElementTransformation(0)), tol_);
}

TEST_F(MFEMHesthaven2D, segmentFromTriangleJacobianO1)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	auto f_trans = m.GetFaceElementTransformations(0);
	ASSERT_EQ(31, f_trans->GetConfigurationMask());
	ASSERT_NEAR(sqrt(2.0), calculateSurfaceJacobian(f_trans->GetFE(), f_trans, fes), tol_);
}

TEST_F(MFEMHesthaven2D, segmentFromTriangleJacobianO2)
{
	const int basis_order = 2;
	auto m{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	auto f_trans = m.GetFaceElementTransformations(0);
	ASSERT_EQ(31, f_trans->GetConfigurationMask());
	ASSERT_NEAR(sqrt(2.0), calculateSurfaceJacobian(f_trans->GetFE(), f_trans, fes), tol_);
}

TEST_F(MFEMHesthaven2D, segmentFromTriangleJacobianO3)
{
	const int basis_order = 3;
	auto m{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	auto f_trans = m.GetFaceElementTransformations(0);
	ASSERT_EQ(31, f_trans->GetConfigurationMask());
	ASSERT_NEAR(sqrt(2.0), calculateSurfaceJacobian(f_trans->GetFE(), f_trans, fes), tol_);
}

TEST_F(MFEMHesthaven2D, connectivityMapO1)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	const auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	GlobalConnectivityMap element_connectivity_map = assembleGlobalConnectivityMap(m, &fec);

	std::vector<std::pair<int, int>> expected_connectivity_pairs({
		{0,4},
		{1,3},
		{1,1},
		{2,2}, 
		{0,0},
		{2,2}, 
		
		{4,0}, 
		{3,1}, 
		{4,4},
		{5,5},
		{3,3},
		{5,5}
	});

	for (auto p{ 0 }; p < expected_connectivity_pairs.size(); p++) {
		EXPECT_EQ(expected_connectivity_pairs[p], element_connectivity_map[p]);
	}
}

TEST_F(MFEMHesthaven2D, connectivityMapO2)
{
	const int basis_order = 2;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	const auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	GlobalConnectivityMap element_connectivity_map = assembleGlobalConnectivityMap(m, &fec);

	std::vector<std::pair<int, int>> expected_connectivity_pairs({ 
		{0,8},   
		{1,7}, 
		{2,6},
		{2,2},
		{4,4},
		{5,5},
		{0,0},
		{3,3},
	    {5,5}, 	
		
		{8,0},
		{7,1},
		{6,2},
		{8,8},
		{10,10},
		{11,11},
		{6,6},
		{9,9},
		{11,11} 
	});

	for (auto p{ 0 }; p < expected_connectivity_pairs.size(); p++) {
		EXPECT_EQ(expected_connectivity_pairs[p], element_connectivity_map[p]);
	}
}

TEST_F(MFEMHesthaven2D, inverseMassMatrixFromSubMeshO1)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };
	auto att_map{ mapOriginalAttributes(m) };
	auto sm{ assembleInteriorFaceSubMesh(m_copy, att_map) };

	// Compute two-element flux matrix 
	auto sm_fes{ FiniteElementSpace(&sm, &fec) };

	auto mass_mat{ assembleMassMatrix(sm_fes) };

	DynamicMatrix el_0_mass_inverse = getElementMassMatrixFromGlobal(0, mass_mat).inverse();
	DynamicMatrix el_1_mass_inverse = getElementMassMatrixFromGlobal(1, mass_mat).inverse();

	DynamicMatrix expectedInverseMass{
		{ 4.5000, -1.5000, -1.5000},
		{-1.5000,  4.5000, -1.5000},
		{-1.5000, -1.5000,  4.5000}
	};

	for (auto r{ 0 }; r < expectedInverseMass.rows(); r++) {
		for (auto c{ 0 }; c < expectedInverseMass.cols(); c++) {
			ASSERT_NEAR(expectedInverseMass(r,c), el_0_mass_inverse(r,c) * hesthaven_triangle_scaling_factor, tol_);
			ASSERT_NEAR(expectedInverseMass(r,c), el_1_mass_inverse(r,c) * hesthaven_triangle_scaling_factor, tol_);
		}
	}
}

TEST_F(MFEMHesthaven2D, inverseMassMatrixFromSubMeshO2)
{
	const int basis_order = 2;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };
	auto att_map{ mapOriginalAttributes(m) };
	auto sm{ assembleInteriorFaceSubMesh(m_copy, att_map) };

	// Compute two-element flux matrix
	auto sm_fes{ FiniteElementSpace(&sm, &fec) };

	auto mass_mat{ assembleMassMatrix(sm_fes) };

	DynamicMatrix el_0_mass_inverse = getElementMassMatrixFromGlobal(0, mass_mat).inverse();
	DynamicMatrix el_1_mass_inverse = getElementMassMatrixFromGlobal(1, mass_mat).inverse();

	DynamicMatrix expectedInverseMass{
		{18.0000, -0.7500,  3.0000, -0.7500,  3.0000,  3.0000},
		{-0.7500,  4.8750, -0.7500, -1.6875, -1.6875,  3.0000},
		{ 3.0000, -0.7500, 18.0000,  3.0000, -0.7500,  3.0000},
		{-0.7500, -1.6875,  3.0000,  4.8750, -1.6875, -0.7500},
		{ 3.0000, -1.6875, -0.7500, -1.6875,  4.8750, -0.7500},
		{ 3.0000,  3.0000,  3.0000, -0.7500, -0.7500, 18.0000}
	};

	for (auto r{ 0 }; r < expectedInverseMass.rows(); r++) {
		for (auto c{ 0 }; c < expectedInverseMass.cols(); c++) {
			ASSERT_NEAR(expectedInverseMass(r, c), el_0_mass_inverse(r, c) * hesthaven_triangle_scaling_factor, tol_);
			ASSERT_NEAR(expectedInverseMass(r, c), el_1_mass_inverse(r, c) * hesthaven_triangle_scaling_factor, tol_);
		}
	}
}

TEST_F(MFEMHesthaven2D, massMatrixTriangleFaceO1)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };

	Array<int> marker;
	marker.Append(hesthavenMeshingTag);
	m_copy.SetAttribute(0, hesthavenMeshingTag);
	auto sm = SubMesh::CreateFromDomain(m_copy, marker);
	FiniteElementSpace sm_fes(&sm, &fec);

	auto boundary_markers = assembleBoundaryMarkers(sm_fes);

	for (auto f{ 0 }; f < sm_fes.GetNF(); f++) {
		sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
	}

	{
		auto surface_matrix{ assembleConnectivityFaceMassMatrix(sm_fes, boundary_markers[0]) };

		DynamicMatrix expected_emat{
			{0.6667, 0.3333},
			{0.3333, 0.6667},
			{0.0000, 0.0000}
		};

		for (auto r{ 0 }; r < surface_matrix.rows(); r++) {
			for (auto c{ 0 }; c < surface_matrix.cols(); c++) {
				EXPECT_NEAR(expected_emat(r, c), surface_matrix(r, c), tol_);
			}
		}

	}

	{
		auto surface_matrix{ assembleConnectivityFaceMassMatrix(sm_fes, boundary_markers[1]) };

		DynamicMatrix expected_emat{
			{0.0000, 0.0000},
			{0.6667, 0.3333},
			{0.3333, 0.6667}
		};

		for (auto r{ 0 }; r < surface_matrix.rows(); r++) {
			for (auto c{ 0 }; c < surface_matrix.cols(); c++) {
				EXPECT_NEAR(expected_emat(r, c), surface_matrix(r, c), tol_);
			}
		}

	}

	{
		auto surface_matrix{ assembleConnectivityFaceMassMatrix(sm_fes, boundary_markers[2]) };

		DynamicMatrix expected_emat{
			{0.6667, 0.3333},
			{0.0000, 0.0000},
			{0.3333, 0.6667}
		};

		for (auto r{ 0 }; r < surface_matrix.rows(); r++) {
			for (auto c{ 0 }; c < surface_matrix.cols(); c++) {
				EXPECT_NEAR(expected_emat(r, c), surface_matrix(r, c), tol_);
			}
		}

	}


}

TEST_F(MFEMHesthaven2D, EmatO1)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };

	Array<int> marker;
	marker.Append(hesthavenMeshingTag);
	m_copy.SetAttribute(0, hesthavenMeshingTag);
	auto sm = SubMesh::CreateFromDomain(m_copy, marker);
	FiniteElementSpace sm_fes(&sm, &fec);

	auto boundary_markers = assembleBoundaryMarkers(sm_fes);

	for (auto f{ 0 }; f < sm_fes.GetNF(); f++) {
		sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
	}

	DynamicMatrix emat = assembleEmat(sm_fes, boundary_markers);

	DynamicMatrix expected_emat{
	   {0.6667, 0.3333, 0.0000, 0.0000, 0.6667, 0.3333},
	   {0.3333, 0.6667, 0.6667, 0.3333, 0.0000, 0.0000},
	   {0.0000, 0.0000, 0.3333, 0.6667, 0.3333, 0.6667}
	};

	for (auto r{ 0 }; r < emat.rows(); r++) {
		for (auto c{ 0 }; c < emat.cols(); c++) {
			EXPECT_NEAR(expected_emat(r, c), emat(r, c), tol_);
		}
	}
}

TEST_F(MFEMHesthaven2D, EmatO2)
{
	const int basis_order = 2;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };

	// Create a singular element submesh.
	Array<int> marker;
	marker.Append(hesthavenMeshingTag);
	m_copy.SetAttribute(0, hesthavenMeshingTag);
	auto sm = SubMesh::CreateFromDomain(m_copy, marker);

	FiniteElementSpace sm_fes(&sm, &fec);
	auto boundary_markers = assembleBoundaryMarkers(sm_fes);

	for (auto f{ 0 }; f < sm_fes.GetNF(); f++) {
		sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
	}

	DynamicMatrix emat = assembleEmat(sm_fes, boundary_markers);

	DynamicMatrix expected_emat{
    { 0.2667, 0.1333, -0.0667,  0.0000, 0.0000,  0.0000,  0.2667, 0.1333, -0.0667},
    { 0.1333, 1.0667,  0.1333,  0.0000, 0.0000,  0.0000,  0.0000, 0.0000,  0.0000},
    {-0.0667, 0.1333,  0.2667,  0.2667, 0.1333, -0.0667,  0.0000, 0.0000,  0.0000},
    { 0.0000, 0.0000,  0.0000,  0.0000, 0.0000,  0.0000,  0.1333, 1.0667,  0.1333},
    { 0.0000, 0.0000,  0.0000,  0.1333, 1.0667,  0.1333,  0.0000, 0.0000,  0.0000},
    { 0.0000, 0.0000,  0.0000, -0.0667, 0.1333,  0.2667, -0.0667, 0.1333,  0.2667}
	};

	std::cout << emat << std::endl;

	for (auto r{ 0 }; r < emat.rows(); r++) {
		for (auto c{ 0 }; c < emat.cols(); c++) {
			EXPECT_NEAR(expected_emat(r, c), emat(r, c), tol_);
		}
	}
}

TEST_F(MFEMHesthaven2D, normals)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };
	auto att_map{ mapOriginalAttributes(m) };

	// Create a singular element submesh
	Array<int> marker;
	marker.Append(hesthavenMeshingTag);

	DynamicMatrix normal_mat_x, normal_mat_y;
	normal_mat_x.resize((basis_order + 1) * m.GetElement(0)->GetNEdges(), fes.GetNE());
	normal_mat_y.resizeLike(normal_mat_x);

	for (auto e{ 0 }; e < fes.GetNE(); e++) {

		m_copy.SetAttribute(e, hesthavenMeshingTag);
		auto sm = SubMesh::CreateFromDomain(m_copy, marker);
		FiniteElementSpace sm_fes(&sm, &fec);
		m_copy.SetAttribute(e, att_map[e]);

		for (auto f{ 0 }; f < sm_fes.GetMesh()->GetNEdges(); f++){
			Vector normal(sm_fes.GetMesh()->SpaceDimension());
			ElementTransformation* f_trans = sm_fes.GetMesh()->GetEdgeTransformation(f);
			f_trans->SetIntPoint(&Geometries.GetCenter(f_trans->GetGeometryType()));
			CalcOrtho(f_trans->Jacobian(), normal);
		
			normal_mat_x( 2 * f     , e) = normal[0] / normal.Norml2();
			normal_mat_x((2 * f) + 1, e) = normal[0] / normal.Norml2();
			normal_mat_y( 2 * f     , e) = normal[1] / normal.Norml2();
			normal_mat_y((2 * f) + 1, e) = normal[1] / normal.Norml2();
		}
	}

	DynamicMatrix expected_normal_x{
	  {-1.0000, -0.7071},
	  {-1.0000, -0.7071},
	  { 0.7071,  0.0000},
	  { 0.7071,  0.0000},
	  { 0.0000,  1.0000},
	  { 0.0000,  1.0000}
	};

	DynamicMatrix expected_normal_y{
	  {-0.0000,  0.7071},
	  {-0.0000,  0.7071},
	  {-0.7071, -1.0000},
	  {-0.7071, -1.0000},
	  { 1.0000,  0.0000},
	  { 1.0000,  0.0000}
	};

	// Face ordering is different between Hesthaven and MFEM
	// Values are manually checked based on ordering

	EXPECT_NEAR(expected_normal_x(0, 0), normal_mat_x(4, 0), tol_);
	EXPECT_NEAR(expected_normal_x(1, 0), normal_mat_x(5, 0), tol_);
	EXPECT_NEAR(expected_normal_x(2, 0), normal_mat_x(0, 0), tol_);
	EXPECT_NEAR(expected_normal_x(3, 0), normal_mat_x(1, 0), tol_);
	EXPECT_NEAR(expected_normal_x(4, 0), normal_mat_x(2, 0), tol_);
	EXPECT_NEAR(expected_normal_x(5, 0), normal_mat_x(3, 0), tol_);

	EXPECT_NEAR(expected_normal_x(0, 1), normal_mat_x(0, 1), tol_);
	EXPECT_NEAR(expected_normal_x(1, 1), normal_mat_x(1, 1), tol_);
	EXPECT_NEAR(expected_normal_x(2, 1), normal_mat_x(2, 1), tol_);
	EXPECT_NEAR(expected_normal_x(3, 1), normal_mat_x(3, 1), tol_);
	EXPECT_NEAR(expected_normal_x(4, 1), normal_mat_x(4, 1), tol_);
	EXPECT_NEAR(expected_normal_x(5, 1), normal_mat_x(5, 1), tol_);

	EXPECT_NEAR(expected_normal_y(0, 0), normal_mat_y(4, 0), tol_);
	EXPECT_NEAR(expected_normal_y(1, 0), normal_mat_y(5, 0), tol_);
	EXPECT_NEAR(expected_normal_y(2, 0), normal_mat_y(0, 0), tol_);
	EXPECT_NEAR(expected_normal_y(3, 0), normal_mat_y(1, 0), tol_);
	EXPECT_NEAR(expected_normal_y(4, 0), normal_mat_y(2, 0), tol_);
	EXPECT_NEAR(expected_normal_y(5, 0), normal_mat_y(3, 0), tol_);

	EXPECT_NEAR(expected_normal_y(0, 1), normal_mat_y(0, 1), tol_);
	EXPECT_NEAR(expected_normal_y(1, 1), normal_mat_y(1, 1), tol_);
	EXPECT_NEAR(expected_normal_y(2, 1), normal_mat_y(2, 1), tol_);
	EXPECT_NEAR(expected_normal_y(3, 1), normal_mat_y(3, 1), tol_);
	EXPECT_NEAR(expected_normal_y(4, 1), normal_mat_y(4, 1), tol_);
	EXPECT_NEAR(expected_normal_y(5, 1), normal_mat_y(5, 1), tol_);

}