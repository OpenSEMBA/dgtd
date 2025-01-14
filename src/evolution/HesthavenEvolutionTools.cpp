#include "evolution/HesthavenEvolutionTools.h"

namespace maxwell {

	using namespace mfem;

	InteriorFaceConnectivityMaps mapConnectivity(const DynamicMatrix& flux_mat)
	{
		// Make map
		std::vector<int> vmapM, vmapP;
		for (auto r{ 0 }; r < flux_mat.rows() / 2; r++) {
			for (auto c{ 0 }; c < flux_mat.cols(); c++) {
				if (r != c && flux_mat(r, c) == flux_mat(r, r) && flux_mat(r, r) > 1e-5) {
					vmapM.emplace_back(r);
					vmapP.emplace_back(c);
				}
			}
		}
		return std::make_pair(vmapM, vmapP);
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

	void restoreOriginalAttributesAfterSubMeshing(FaceElementTransformations* f_trans, Mesh& m_copy, const std::map<int, Attribute>& att_map)
	{
		m_copy.SetAttribute(f_trans->Elem1No, att_map.at(f_trans->Elem1No));
		if (f_trans->Elem2No) {
			m_copy.SetAttribute(f_trans->Elem2No, att_map.at(f_trans->Elem2No));
		}
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
		m_copy.SetAttribute(f_trans->Elem1No, hesthavenMeshingTag);
		m_copy.SetAttribute(f_trans->Elem2No, hesthavenMeshingTag);
	}

	DynamicMatrix loadMatrixWithValues(const DynamicMatrix& global, const int start_row, const int start_col)
	{
		DynamicMatrix res;
		res.resize(int(global.rows() / 2), int(global.cols() / 2));
		for (auto r{ start_row }; r < start_row + int(global.rows() / 2); r++) {
			for (auto c{ start_col }; c < start_col + int(global.cols() / 2); c++) {
				res(r - start_row, c - start_col) = global(r, c);
			}
		}
		return res;
	}

	DynamicMatrix getElementMassMatrixFromGlobal(const int el, const DynamicMatrix& global)
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

	DynamicMatrix getFaceMassMatrixFromGlobal(const DynamicMatrix& global)
	{
		DynamicMatrix res;
		res.resize(int(global.rows() / 2) - 1, int(global.cols() / 2) - 1);
		for (auto r{ 0 }; r < int(global.rows() / 2) - 1; r++) {
			for (auto c{ 0 }; c < int(global.cols() / 2) - 1; c++) {
				res(r, c) = global(r, c);
			}
		}
		return res;
	}

	void removeColumn(DynamicMatrix& matrix, const int colToRemove)
	{
		auto numRows = matrix.rows();
		auto numCols = matrix.cols() - 1;

		if (colToRemove < numCols)
			matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

		matrix.conservativeResize(numRows, numCols);
	}

	void removeZeroColumns(DynamicMatrix& matrix)
	{
		std::vector<int> cols_to_del;
		for (auto c{ 0 }; c < matrix.cols(); c++) {
			bool isZero = true;
			for (auto r{ 0 }; r < matrix.rows(); r++) {
				if (matrix(r, c) != 0.0)
				{
					isZero = false;
					break;
				}
			}
			if (isZero == true) {
				cols_to_del.emplace_back(c);
			}
		}
		for (auto it = cols_to_del.rbegin(); it != cols_to_del.rend(); ++it) {
			removeColumn(matrix, *it);
		}
	}

	std::vector<Array<int>> assembleBoundaryMarkers(FiniteElementSpace& fes)
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

	DynamicMatrix assembleConnectivityFaceMassMatrix(FiniteElementSpace& fes, Array<int> boundary_marker)
	{
		auto bf = assembleFaceMassBilinearForm(fes, boundary_marker);
		auto res = toEigen(*bf->SpMat().ToDenseMatrix());
		removeZeroColumns(res);
		return res;
	}

	DynamicMatrix assembleEmat(FiniteElementSpace& fes, std::vector<Array<int>>& boundary_markers)
	{
		DynamicMatrix res;
		res.resize(fes.GetNDofs(), fes.GetNF() * (fes.GetMaxElementOrder() + 1));
		for (auto f{ 0 }; f < fes.GetNF(); f++) {
			auto surface_matrix{ assembleConnectivityFaceMassMatrix(fes, boundary_markers[f]) };
			res.block(0, f * (fes.GetMaxElementOrder() + 1), fes.GetNDofs(), fes.GetMaxElementOrder() + 1) = surface_matrix;
		}
		return res;
	}

	SubMesh assembleInteriorFaceSubMesh(Mesh& m, const FaceElementTransformations& trans)
	{
		m.SetAttribute(trans.Elem1No, hesthavenMeshingTag);
		m.SetAttribute(trans.Elem2No, hesthavenMeshingTag);
		Array<int> volume_marker;
		volume_marker.Append(hesthavenMeshingTag);
		return SubMesh::CreateFromDomain(m, volume_marker);
	}

	DynamicMatrix assembleInteriorFluxMatrix(FiniteElementSpace& fes)
	{
		BilinearForm bf(&fes);
		bf.AddInteriorFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0));
		bf.Assemble();
		bf.Finalize();
		return toEigen(*bf.SpMat().ToDenseMatrix());
	}

	DynamicMatrix assembleBoundaryFluxMatrix(FiniteElementSpace& fes)
	{
		BilinearForm bf(&fes);
		bf.AddBdrFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0));
		bf.Assemble();
		bf.Finalize();
		return toEigen(*bf.SpMat().ToDenseMatrix());
	}

	void appendConnectivityMapsFromInteriorFace(const FaceElementTransformations& trans, const int element_index, FiniteElementSpace& fes, GlobalConnectivityMap& map)
	{
		auto int_flux_mat{ assembleInteriorFluxMatrix(fes) };
		auto maps{ mapConnectivity(int_flux_mat) };
		for (auto index{ 0 }; index < maps.first.size(); index++) {
			if (trans.Elem1No == element_index) {
				map.push_back(std::pair(maps.first[index], maps.second[index]));
			}
			else if (trans.Elem2No == element_index) {
				map.push_back(std::pair(maps.second[index], maps.first[index]));
			}
			else {
				throw std::runtime_error("Wrong element number in connectivity map assembly.");
			}
		}
	}

	void appendConnectivityMapsFromBoundaryFace(FiniteElementSpace& fes, FiniteElementSpace& sm_fes, const DynamicMatrix& surface_matrix, GlobalConnectivityMap& map)
	{
		GridFunction gf_parent(&fes);
		GridFunction gf_son(&sm_fes);
		auto mesh = dynamic_cast<SubMesh*>(sm_fes.GetMesh());
		auto transfer_map = mesh->CreateTransferMap(gf_son, gf_parent);
		ConstantCoefficient zero(0.0);
		for (auto r{ 0 }; r < surface_matrix.rows(); r++) {
			if (surface_matrix(r, 0) != 0.0) {
				gf_son.ProjectCoefficient(zero);
				gf_parent.ProjectCoefficient(zero);
				gf_son[r] = 1.0;
				transfer_map.Transfer(gf_son, gf_parent);
				for (auto v{ 0 }; v < gf_parent.Size(); v++) {
					if (gf_parent[v] != 0.0) {
						map.push_back(std::make_pair(v, v));
						break;
					}
				}
			}
		}
	}

	void tagBdrAttributesForSubMesh(const int edge, SubMesh& sm)
	{
		for (auto b{ 0 }; b < sm.bdr_attributes.Size(); b++) {
			sm.SetBdrAttribute(b, 1);
		}
		sm.SetBdrAttribute(edge, hesthavenMeshingTag);
		sm.bdr_attributes.Append(hesthavenMeshingTag);
	}

	GlobalConnectivityMap assembleGlobalConnectivityMap(Mesh& m, const L2_FECollection* fec)
	{
		GlobalConnectivityMap res;
		FiniteElementSpace fes(&m, fec);

		std::map<FaceId, bool> global_face_is_interior;
		for (auto f{ 0 }; f < m.GetNEdges(); f++) {
			global_face_is_interior[f] = m.FaceIsInterior(f);
		}

		Table global_element_to_edge = m.ElementToEdgeTable();

		Array<int> volume_marker;
		volume_marker.Append(hesthavenMeshingTag);

		Array<int> boundary_marker(hesthavenMeshingTag);
		boundary_marker = 0;
		boundary_marker[hesthavenMeshingTag - 1] = 1;

		for (auto e{ 0 }; e < m.GetNE(); e++) {

			Array<int> local_edge_index_to_global_edge_index;
			global_element_to_edge.GetRow(e, local_edge_index_to_global_edge_index);

			for (auto local_edge{ 0 }; local_edge < m.GetElement(e)->GetNEdges(); local_edge++) {

				auto m_copy{ Mesh(m) };

				if (!global_face_is_interior[local_edge_index_to_global_edge_index[local_edge]]) {

					m_copy.SetAttribute(e, hesthavenMeshingTag);
					auto sm = SubMesh::CreateFromDomain(m_copy, volume_marker);
					tagBdrAttributesForSubMesh(local_edge, sm);
					FiniteElementSpace sm_fes(&sm, fec);
					appendConnectivityMapsFromBoundaryFace(
						fes,
						sm_fes,
						assembleConnectivityFaceMassMatrix(sm_fes, boundary_marker),
						res);

				}
				else {

					FaceElementTransformations* f_trans = m_copy.GetFaceElementTransformations(local_edge_index_to_global_edge_index[local_edge]);
					auto sm = assembleInteriorFaceSubMesh(m_copy, *f_trans);
					FiniteElementSpace sm_fes(&sm, fec);
					appendConnectivityMapsFromInteriorFace(
						*f_trans,
						e,
						sm_fes,
						res);

				}
			}
		}
		return res;
	}
}