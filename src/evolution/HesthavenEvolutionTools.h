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
	using BdrCondToNodes = std::pair<BdrCond, std::vector<int>>;
	using GlobalBoundaryMap = std::vector<BdrCondToNodes>;
	using ConnectivityMap = std::vector<std::pair<int, int>>;
	using Emat = std::vector<const DynamicMatrix*>;
	using Directional = std::array<const DynamicMatrix*, 3>;
	using Normals = std::array<Eigen::VectorXd, 3>;

	struct MatrixCompareLessThan {

		bool operator()(const Eigen::MatrixXd& lhs, const Eigen::MatrixXd& rhs) const {

			for (int i{ 0 }; i < lhs.rows(); i++) {
				for (int j{ 0 }; j < lhs.cols(); j++) {
					if (lhs(i, j) != rhs(i, j)) {
						return lhs(i, j) < rhs(i, j);
					}
				}
			}
			return false;
		}

	};

	using StorageIterator = std::set<DynamicMatrix, MatrixCompareLessThan>::iterator;

	struct HesthavenElement {
		ElementId id;
		Geometry::Type geom;
		Directional dir;
		Emat emat;
		const DynamicMatrix* invmass;
		Normals normals;
		Eigen::VectorXd fscale;
	};

	struct BoundaryElement {
		Array<FaceId> globalFaceIndex;
		Array<Attribute> faceAttributes;
		Array<ElementId> elementId;
	};

	struct ConnectivityMaps {
		std::vector<int> vmapM;
		std::vector<int> vmapP;
		std::vector<int> vmapB;
	};

	void restoreOriginalAttributesAfterSubMeshing(FaceElementTransformations* faceTrans, Mesh&, const std::map<int, Attribute>&);
	void restoreOriginalAttributesAfterSubMeshing(ElementId, Mesh&, const std::map<int, Attribute>&);
	void markElementsForSubMeshing(FaceElementTransformations* faceTrans, Mesh&);
	void removeColumn(DynamicMatrix&, const int colToRemove);
	void removeZeroColumns(DynamicMatrix&);
	void appendConnectivityMapsFromInteriorFace(const FaceElementTransformations&, const ElementId, FiniteElementSpace&, GlobalConnectivityMap&);
	void appendConnectivityMapsFromBoundaryFace(FiniteElementSpace& globalFES, FiniteElementSpace& submeshFES, const DynamicMatrix& surfaceMatrix, GlobalConnectivityMap&);
	void tagBdrAttributesForSubMesh(const FaceId, SubMesh& sm);

	std::map<int, Attribute> mapOriginalAttributes(const Mesh& m);

	DynamicMatrix loadMatrixWithValues(const DynamicMatrix& global, const int startRow, const int startCol);
	DynamicMatrix getElementMassMatrixFromGlobal(const int el, const DynamicMatrix& global);
	DynamicMatrix getFaceMassMatrixFromGlobal(const DynamicMatrix& global);
	DynamicMatrix assembleConnectivityFaceMassMatrix(FiniteElementSpace& fes, Array<int> boundaryMarker);
	DynamicMatrix assembleEmat(FiniteElementSpace& fes, std::vector<Array<int>>& boundaryMarkers);
	DynamicMatrix assembleInteriorFluxMatrix(FiniteElementSpace&);
	DynamicMatrix assembleBoundaryFluxMatrix(FiniteElementSpace&);
	
	SubMesh assembleInteriorFaceSubMesh(Mesh& m, const FaceElementTransformations& trans);
	
	std::vector<Array<int>> assembleBoundaryMarkers(FiniteElementSpace&);
	std::unique_ptr<BilinearForm> assembleFaceMassBilinearForm(FiniteElementSpace&, Array<int>& boundaryMarker);
	InteriorFaceConnectivityMaps mapConnectivity(const DynamicMatrix& fluxMatrix);
	Array<int> getFacesForElement(const Mesh&, const ElementId);
	FaceElementTransformations* getInteriorFaceTransformation(Mesh&, const Array<int>& faces);
	GlobalConnectivityMap assembleGlobalConnectivityMap(Mesh, const L2_FECollection*);



}