#pragma once

#include <unordered_set>

#include "mfemExtension/BilinearIntegrators.h"
#include "math/EigenMfemTools.h"
#include "components/Types.h"

namespace maxwell {

	using namespace mfem;

	using InteriorFaceConnectivityMaps = std::pair<std::vector<int>, std::vector<int>>;
	using DynamicMatrix = Eigen::MatrixXd;
	using GlobalConnectivityMap = std::vector<std::pair<int, int>>;
	using ConnectivityMap = std::vector<std::pair<int, int>>;

	struct Normals {
		std::vector<double> X;
		std::vector<double> Y;
		std::vector<double> Z;
	};

	struct Emat {
		const DynamicMatrix* face0;
		const DynamicMatrix* face1;
		const DynamicMatrix* face2;
	};

	struct HesthavenElement {
		ElementId id;
		Geometry::Type geom;
		Emat emat;
		Normals normals;
		std::vector<double> fscale;
		ConnectivityMap connect;
	};

	struct ConnectivityMaps {
		std::vector<int> vmapM;
		std::vector<int> vmapP;
		std::vector<int> vmapB;
	};

	enum FacesPerGeom {
		TRIANGLE = 3,
		QUADRILATERAL = 4,
		TETRAHEDRON = 4,
		HEXAHEDRON = 6
	};

	InteriorFaceConnectivityMaps mapConnectivity(const DynamicMatrix& flux_mat);
	void restoreOriginalAttributesAfterSubMeshing(FaceElementTransformations* f_trans, Mesh& m_copy, const std::map<int, Attribute>& att_map);
	Array<int> getFacesForElement(const Mesh& m_copy, const int el);
	FaceElementTransformations* getInteriorFaceTransformation(Mesh& m_copy, const Array<int>& faces);
	void markElementsForSubMeshing(FaceElementTransformations* f_trans, Mesh& m_copy);
	SubMesh createSubMeshFromInteriorFace(Mesh& m_copy, const std::map<int, Attribute>& att_map);
	DynamicMatrix loadMatrixWithValues(const DynamicMatrix& global, const int start_row, const int start_col);
	DynamicMatrix getElementMassMatrixFromGlobal(const int el, const DynamicMatrix& global);
	DynamicMatrix getFaceMassMatrixFromGlobal(const DynamicMatrix& global);
	void removeColumn(DynamicMatrix& matrix, const int colToRemove);
	void removeZeroColumns(DynamicMatrix& matrix);
	std::vector<Array<int>> assembleBoundaryMarkers(FiniteElementSpace& fes);
	std::unique_ptr<BilinearForm> assembleFaceMassBilinearForm(FiniteElementSpace& fes, Array<int>& boundary_marker);
	DynamicMatrix assembleConnectivityFaceMassMatrix(FiniteElementSpace& fes, Array<int> boundary_marker);
	DynamicMatrix assembleEmat(FiniteElementSpace& fes, std::vector<Array<int>>& boundary_markers);
	SubMesh assembleInteriorFaceSubMesh(Mesh& m, const FaceElementTransformations& trans);
	DynamicMatrix assembleInteriorFluxMatrix(FiniteElementSpace& fes);
	DynamicMatrix assembleBoundaryFluxMatrix(FiniteElementSpace& fes);
	void appendConnectivityMapsFromInteriorFace(const FaceElementTransformations& trans, const int element_index, FiniteElementSpace& fes, GlobalConnectivityMap& map);
	void appendConnectivityMapsFromBoundaryFace(FiniteElementSpace& fes, FiniteElementSpace& sm_fes, const DynamicMatrix& surface_matrix, GlobalConnectivityMap& map);
	void tagBdrAttributesForSubMesh(const int edge, SubMesh& sm);
	GlobalConnectivityMap assembleGlobalConnectivityMap(Mesh& m, const L2_FECollection* fec);



}