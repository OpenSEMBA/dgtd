#pragma once

#include "mfemExtension/BilinearIntegrators.h"
#include "math/EigenMfemTools.h"

using namespace mfem;

namespace maxwell {

	Attribute hesthaven_element_tag_ = 777;

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

	SubMesh createSubMeshFromInteriorFace(Mesh& m_copy, const std::map<int, Attribute>& att_map)
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
		for (auto r{ start_row }; r < start_row + int(global.rows() / 2); r++) {
			for (auto c{ start_col }; c < start_col + int(global.cols() / 2); c++) {
				res(r - start_row, c - start_col) = global(r, c);
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

	std::unique_ptr<BilinearForm> assembleFaceMassBilinearForm(FiniteElementSpace& fes, Array<int>& boundary_marker)
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
}