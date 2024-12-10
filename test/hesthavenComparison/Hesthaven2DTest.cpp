#include <gtest/gtest.h>

#include "TestUtils.h"
#include "HesthavenFunctions.h"
#include "math/EigenMfemTools.h"
#include "evolution/EvolutionMethods.h"
#include "mfemExtension/BilinearIntegrators.h"

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
	Attribute hesthaven_element_tag_ = 777;
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

	std::map<int, Attribute> mapOriginalAttributes(const Mesh& m)
	{
		// Create backup of attributes
		auto res{ std::map<int,Attribute>() };
		for (auto e{ 0 }; e < m.GetNE(); e++) {
			res[e] = m.GetAttribute(e);
		}
		return res;
	}
	Eigen::MatrixXd assembleMassMatrix(FiniteElementSpace& fes)
	{
		BilinearForm bf(&fes);
		ConstantCoefficient one(1.0);
		bf.AddDomainIntegrator(new MassIntegrator(one));
		bf.Assemble();
		bf.Finalize();

		return toEigen(*bf.SpMat().ToDenseMatrix());
	}


	Eigen::MatrixXd assembleFluxMatrix(FiniteElementSpace& fes)
	{
		BilinearForm bf(&fes);
		bf.AddInteriorFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0));
		bf.Assemble();
		bf.Finalize();
		return toEigen(*bf.SpMat().ToDenseMatrix());
	}

	std::map<int, int> mapConnectivityLocal(const Eigen::MatrixXd& flux_mat)
	{
		// Make map
		std::map<int, int> res;
		for (auto r{ 0 }; r < flux_mat.rows(); r++) {
			for (auto c{ r + 1 }; c < flux_mat.cols(); c++) {
				if (flux_mat(r, c) == flux_mat(r, r) && flux_mat(r, r) > 1e-5) {
					res[r] = c;
				}
			}
		}
		return res;
	}

	void restoreOriginalAttributesAfterSubMeshing(const FaceElementTransformations* f_trans, Mesh& m_copy, const std::map<int, Attribute>& att_map)
	{
		m_copy.SetAttribute(f_trans->Elem1No, att_map.at(f_trans->Elem1No));
		m_copy.SetAttribute(f_trans->Elem2No, att_map.at(f_trans->Elem2No));
	}

	Array<int> getFacesForElement(const Mesh& m_copy, const int el)
	{
		auto e2f = m_copy.ElementToEdgeTable();
		e2f.Finalize();
		Array<int> res;
		e2f.GetRow(el, res);
		return res;
	}

	FaceElementTransformations* getInteriorFaceTransformation(Mesh& m_copy, const Array<int>& faces)
	{
		for (auto f{ 0 }; f < faces.Size(); f++) {
			FaceElementTransformations* res = m_copy.GetFaceElementTransformations(faces[f]);
			if (res->GetConfigurationMask() == 31) {
				return res;
			}
		}
		throw std::runtime_error("There is no interior face for the selected element in getInteriorFaceTransformation.");
	}

	void markElementsForSubMeshing(FaceElementTransformations* f_trans, Mesh& m_copy)
	{
		m_copy.SetAttribute(f_trans->Elem1No, hesthaven_element_tag_);
		m_copy.SetAttribute(f_trans->Elem2No, hesthaven_element_tag_);
	}

	SubMesh createSubMeshFromInteriorFace(Mesh& m_copy, const std::map<int,Attribute>& att_map)
	{
		Array<int> sm_tag;
		sm_tag.Append(hesthaven_element_tag_);
		auto faces = getFacesForElement(m_copy, 0);
		auto f_trans = getInteriorFaceTransformation(m_copy, faces);
		markElementsForSubMeshing(f_trans, m_copy);
		auto res = SubMesh::CreateFromDomain(m_copy, sm_tag);
		restoreOriginalAttributesAfterSubMeshing(f_trans, m_copy, att_map);
		return res;
	}

	Eigen::MatrixXd loadMatrixWithValues(const Eigen::MatrixXd& global, const int start_row, const int start_col)
	{
		Eigen::MatrixXd res;
		res.resize(int(global.rows() / 2), int(global.cols() / 2));
		for (auto r{ start_row }; r < start_row + int(global.rows()/2); r++){
			for (auto c{ start_col }; c < start_col + int(global.cols() / 2); c++) {
				res(r-start_row,c-start_col) = global(r, c);
			}
		}
		return res;
	}

	Eigen::MatrixXd getElementMassMatrixFromGlobal(const int el, const Eigen::MatrixXd& global)
	{
		switch (el) {
		case 0:
			return loadMatrixWithValues(global, 0, 0);
		case 1:
			return loadMatrixWithValues(global, int(global.rows() / 2), int(global.cols() / 2));
		default:
			throw std::runtime_error("Incorrect element index for getElementMassMatrixFromGlobal");
		}
	}

	Eigen::MatrixXd getFaceMassMatrixFromGlobal(const Eigen::MatrixXd& global)
	{
		Eigen::MatrixXd res;
		res.resize(int(global.rows() / 2) - 1, int(global.cols() / 2) - 1);
		for (auto r{ 0 }; r < int(global.rows() / 2) - 1; r++) {
			for (auto c{ 0 }; c < int(global.cols() / 2) - 1; c++) {
				res(r, c) = global(r, c);
			}
		}
		return res;
	}

	Eigen::MatrixXd assembleEmat(const FiniteElementSpace& fes, const Eigen::MatrixXd& faceMat)
	{
		Eigen::MatrixXd res;
		res.resize(int(fes.GetNDofs()/fes.GetNE()), int(fes.GetNFDofs() / (fes.GetNF() * fes.GetNE())) * int(fes.GetNF() / fes.GetNE()));


	}

	void removeColumn(Eigen::MatrixXd& matrix, const int colToRemove)
	{
		auto numRows = matrix.rows();
		auto numCols = matrix.cols() - 1;

		if (colToRemove < numCols)
			matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

		matrix.conservativeResize(numRows, numCols);
	}

	void removeZeroColumns(Eigen::MatrixXd& matrix)
	{
		for (auto c{ 0 }; c < matrix.cols(); c++) {
			bool isZero = true;
			for (auto r{ 0 }; r < matrix.rows(); r++) {
				if (matrix(r, c) != 0.0)
				{
					isZero = false;
				}
			}
			if (isZero == true) {
				removeColumn(matrix, c);
			}
		}
	}

	std::vector<Array<int>> assembleBoundaryMarkers(const FiniteElementSpace& fes) 
	{		
		std::vector<Array<int>> res;
		for (auto f{ 0 }; f < fes.GetNF(); f++) {
			Array<int> bdr_marker;
			bdr_marker.SetSize(fes.GetMesh()->bdr_attributes.Max());
			bdr_marker = 0;
			bdr_marker[fes.GetMesh()->bdr_attributes[f] - 1] = 1;
			res.emplace_back(bdr_marker);
		}
		return res;
	}

	FiniteElementOperator assembleFaceMassBilinearForm(FiniteElementSpace& fes, Array<int>& boundary_marker) 
	{
		auto res = std::make_unique<BilinearForm>(&fes);
		res->AddBdrFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0), boundary_marker);
		res->Assemble();
		res->Finalize();

		return res;
	}

	Eigen::MatrixXd assembleEmat(FiniteElementSpace& fes, std::vector<Array<int>>& boundary_markers)
	{
		Eigen::MatrixXd res;
		res.resize(fes.GetNDofs(), fes.GetNF() * (fes.GetMaxElementOrder() + 1));
		for (auto f{ 0 }; f < fes.GetNF(); f++) {
			auto bf = assembleFaceMassBilinearForm(fes, boundary_markers[f]);
			auto surface_matrix = toEigen(*bf->SpMat().ToDenseMatrix());
			removeZeroColumns(surface_matrix);
			res.block(0, f * (fes.GetMaxElementOrder() + 1), fes.GetNDofs(), fes.GetMaxElementOrder() + 1) = surface_matrix;
		}
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
	Model model(mesh,GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()));

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

TEST_F(MFEMHesthaven2D, connectivityMap)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };
	auto att_map{ mapOriginalAttributes(m) };
	auto sm{ createSubMeshFromInteriorFace(m_copy, att_map) };

	// Compute two-element flux matrix
	auto sm_fes{ FiniteElementSpace(&sm, &fec) };
	auto flux_mat{ assembleFluxMatrix(sm_fes) };

	// Make connectivity map
	auto map{ mapConnectivityLocal(flux_mat) };
	ASSERT_EQ(4, map.at(0));
	ASSERT_EQ(3, map.at(1));
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
	auto sm{ createSubMeshFromInteriorFace(m_copy, att_map) };

	// Compute two-element flux matrix 
	auto sm_fes{ FiniteElementSpace(&sm, &fec) };

	auto mass_mat{ assembleMassMatrix(sm_fes) };

	Eigen::MatrixXd el_0_mass_inverse = getElementMassMatrixFromGlobal(0, mass_mat).inverse();
	Eigen::MatrixXd el_1_mass_inverse = getElementMassMatrixFromGlobal(1, mass_mat).inverse();

	Eigen::Matrix3d expectedInverseMass{
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
	auto sm{ createSubMeshFromInteriorFace(m_copy, att_map) };

	// Compute two-element flux matrix
	auto sm_fes{ FiniteElementSpace(&sm, &fec) };

	auto mass_mat{ assembleMassMatrix(sm_fes) };

	Eigen::MatrixXd el_0_mass_inverse = getElementMassMatrixFromGlobal(0, mass_mat).inverse();
	Eigen::MatrixXd el_1_mass_inverse = getElementMassMatrixFromGlobal(1, mass_mat).inverse();

	Eigen::MatrixXd expectedInverseMass{
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

TEST_F(MFEMHesthaven2D, Emat)
{
	const int basis_order = 1;
	auto m{ Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE) };
	auto fec{ L2_FECollection(basis_order, 2, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	// Not touching the original mesh.
	auto m_copy{ Mesh(m) };
	auto att_map{ mapOriginalAttributes(m) };

	Array<int> marker;
	marker.Append(hesthaven_element_tag_);
	m_copy.SetAttribute(0, hesthaven_element_tag_);
	auto sm = SubMesh::CreateFromDomain(m_copy, marker);
	for (auto f{ 0 }; f < sm.GetElement(0)->GetNEdges(); f++) {
		sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
	}

	FiniteElementSpace sm_fes(&sm, &fec);
	auto boundary_markers = assembleBoundaryMarkers(sm_fes);
	auto emat = assembleEmat(sm_fes, boundary_markers);

	Eigen::MatrixXd expected_emat{
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