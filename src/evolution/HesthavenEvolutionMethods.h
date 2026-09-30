#pragma once

#include <cmath>
#include <unordered_set>

#include "mfemExtension/BilinearIntegrators.h"
#include "solver/SourcesManager.h"
#include "math/EigenMfemTools.h"
#include "components/Model.h"

namespace maxwell {

	using namespace mfem;

	using NodeId = int;
	using Volume = double;
	using Nodes = std::vector<NodeId>;
	using DynamicMatrix = Eigen::MatrixXd;
	using InteriorFaceConnectivityMaps = std::pair<Nodes, Nodes>;
	using BoundaryFaceConnectivityMaps = InteriorFaceConnectivityMaps;
	using NodePair = std::pair<NodeId, NodeId>;
	using ConnectivityVector = std::vector<NodePair>;
	using GlobalConnectivity = ConnectivityVector;
	using BdrCondToNodes = std::pair<BdrCond, Nodes>;
	using Emat = std::vector<const DynamicMatrix*>;
	using Directional = std::array<const DynamicMatrix*, 3>;
	using Normals = std::array<Eigen::VectorXd, 3>;

	struct MatrixCompareLessThan {
		static constexpr double entryTolerance = 1e-10;

		static long long quantize(const double value) {
			return std::llround(value / entryTolerance);
		}

		bool operator()(const Eigen::MatrixXd& lhs, const Eigen::MatrixXd& rhs) const  {
			if (lhs.rows() < rhs.rows()) {
				return true;
			}
			else if (lhs.rows() > rhs.rows()) {
				return false;
			}

			if (lhs.cols() < rhs.cols()) {
				return true;
			}
			else if (lhs.cols() > rhs.cols()) {
				return false;
			}

			for (int i{ 0 }; i < lhs.rows(); i++) {
				for (int j{ 0 }; j < lhs.cols(); j++) {
					const long long lhs_q = quantize(lhs(i, j));
					const long long rhs_q = quantize(rhs(i, j));
					if (lhs_q != rhs_q) {
						return lhs_q < rhs_q;
					}
				}
			}
			return false;
		}

	};

	using StorageIterator = std::set<DynamicMatrix, MatrixCompareLessThan>::iterator;

	struct HesthavenCurvedElement {
		ElementId id;
		Element::Type type;
		Array<int> dofs;
		SparseMatrix matrix;
	};

	struct HesthavenElement {
		ElementId id;
		Element::Type type;
		Volume vol;
		Directional dir;
		Normals normals;
		Eigen::VectorXd fscale;
	};

	struct TrueBoundaryMaps {
		std::vector<Nodes> mapB;
		std::vector<Nodes> vmapB;
	};

	struct InteriorBoundaryMaps {
		std::vector<Nodes> mapBElem1;
		std::vector<Nodes> mapBElem2;
		std::vector<Nodes> vmapBElem1;
		std::vector<Nodes> vmapBElem2;
	};

	struct TotalFieldScatteredFieldMaps {
		std::vector<Nodes> mapBSF;
		std::vector<Nodes> mapBTF;
		std::vector<Nodes> vmapBSF;
		std::vector<Nodes> vmapBTF;
	};

	struct BoundaryMaps {
		TrueBoundaryMaps PEC;
		TrueBoundaryMaps PMC;
		TrueBoundaryMaps SMA;
		InteriorBoundaryMaps intPEC;
		InteriorBoundaryMaps intPMC;
		InteriorBoundaryMaps intSMA;
		TotalFieldScatteredFieldMaps TFSF;
	};

	/// Per global jump slot: 1 if connectivity_.global[v].second is a FaceNbrData index.
	using GlobalPlusIsFaceNbr = std::vector<uint8_t>;

	class Connectivities {
	public:

		Connectivities(Model& model, const FiniteElementSpace& fes);

		GlobalConnectivity global;
		GlobalPlusIsFaceNbr global_plus_is_face_nbr;
		BoundaryMaps boundary;

	private:

		void initBdrConnectivityMaps(const std::vector<Nodes>& bdr2nodes, const std::map<bool, std::vector<BdrElementId>>& isInteriorMap);
		void initIntFaceConnectivityMaps(const BoundaryToMarker& markers);
		void loadIntBdrConditions(const InteriorFaceConnectivityMaps& mapB, const InteriorFaceConnectivityMaps&, const BdrCond&, const double faceOri);
		InteriorFaceConnectivityMaps initInteriorFacesMapB(const InteriorFaceConnectivityMaps& nodePairs) const;

		Model& model_;
		const FiniteElementSpace& fes_;

	};

	struct BoundaryElement {
		Array<FaceId> globalFaceIndex;
		Array<Attribute> faceAttributes;
		Array<ElementId> elementId;
	};

	struct HesthavenFields {
		HesthavenFields(int size);
		std::array<Eigen::VectorXd, 3> e_, h_;
	};

	struct FieldsElementMaps {
		FieldsElementMaps(const Vector& in, FiniteElementSpace&, const ElementId& id);
		std::vector<Eigen::VectorXd> e_, h_;
	};

	struct FieldsInputMaps {
		FieldsInputMaps(const Vector& in, FiniteElementSpace&);
		std::vector<Eigen::Map<Eigen::VectorXd>> e_, h_;
	};	

	struct FieldsOutputMaps {
		FieldsOutputMaps(Vector& out, FiniteElementSpace&);
		std::vector<Eigen::Map<Eigen::VectorXd>> e_, h_;
	};

	struct HesthavenElementJumps {
		HesthavenElementJumps(HesthavenFields& in, ElementId id, Eigen::Index& elFluxSize);
		std::vector<Eigen::Map<Eigen::VectorXd>> e_, h_;
	};

	void restoreOriginalAttributesAfterSubMeshing(FaceElementTransformations* faceTrans, Mesh&, const std::map<int, Attribute>&);
	void restoreOriginalAttributesAfterSubMeshing(ElementId, Mesh&, const std::map<int, Attribute>&);
	void markElementsForSubMeshing(FaceElementTransformations* faceTrans, Mesh&);
	void removeColumn(DynamicMatrix&, const int colToRemove);
	void removeZeroColumns(DynamicMatrix&);
	void appendConnectivityMapsForInteriorFace(const FaceElementTransformations&, FiniteElementSpace& globalFES, FiniteElementSpace& smFES, GlobalConnectivity&, ElementId);
	void appendConnectivityMapsForBoundaryFace(FiniteElementSpace& globalFES, FiniteElementSpace& smFES, const DynamicMatrix& surfaceMatrix, GlobalConnectivity&);
	void tagBdrAttributesForSubMesh(const FaceId, SubMesh& sm);

	const int getNumFaces(const Geometry::Type&); 
	const int getFaceNodeNumByGeomType(const FiniteElementSpace& fes);

	std::map<int, Attribute> mapOriginalAttributes(const Mesh& m);

	DynamicMatrix loadMatrixForTris(const DynamicMatrix& global, const int startRow, const int startCol);
	DynamicMatrix getElemMassMatrixFromGlobal(const int el, const DynamicMatrix& global, const Element::Type elType);
	DynamicMatrix getFaceMassMatrixFromGlobal(const DynamicMatrix& global);
	DynamicMatrix assembleConnectivityFaceMassMatrix(FiniteElementSpace&, Array<int> boundaryMarker);
	DynamicMatrix assembleEmat(FiniteElementSpace&, std::vector<Array<int>>& boundaryMarkers);
	std::unique_ptr<BilinearForm> assembleInteriorFluxMatrix(FiniteElementSpace&);
	DynamicMatrix assembleBoundaryFluxMatrix(FiniteElementSpace&);
	
	const std::vector<Source::Position> buildDoFPositions(const FiniteElementSpace&);

	SubMesh assembleBoundaryFaceSubMesh(Mesh& m, const FaceElementTransformations& trans, const FaceToGeomTag& attMap);
	SubMesh assembleInteriorFaceSubMesh(Mesh& m, const FaceElementTransformations& trans, const FaceToGeomTag& attMap);
	Mesh buildHesthavenRefTetrahedra();

	const std::map<bool, std::vector<BdrElementId>> assembleInteriorOrTrueBdrMap(const FiniteElementSpace& fes);
	std::vector<Nodes> assembleNodeVectorPerBdrFace(std::vector<Array<int>>& bdrNodeMarkers, FiniteElementSpace& , const std::map<bool, std::vector<BdrElementId>>& isInteriorMap);

	BoundaryFaceConnectivityMaps buildConnectivityForBdrFace(const FaceElementTransformations& trans, const FiniteElementSpace& globalFES, FiniteElementSpace& smFES);
	InteriorFaceConnectivityMaps buildConnectivityForInteriorBdrFace(const FaceElementTransformations&, const FiniteElementSpace& globalFES, FiniteElementSpace& smFES);
	std::vector<Array<int>> assembleBoundaryMarkers(FiniteElementSpace&);
	std::unique_ptr<BilinearForm> assembleFaceMassBilinearForm(FiniteElementSpace&, Array<int>& boundaryMarker);
	InteriorFaceConnectivityMaps mapConnectivity(const BilinearForm* fluxMatrix);
	Array<int> getFacesForElement(const Mesh&, const ElementId);
	FaceElementTransformations* getInteriorFaceTransformation(Mesh&, const Array<int>& faces);
	GlobalConnectivity assembleGlobalConnectivityMap(FiniteElementSpace& globalFES,
		const L2_FECollection* fec,
		GlobalPlusIsFaceNbr* plus_is_nbr = nullptr);
	void appendSharedFaceConnectivityForElement(ParFiniteElementSpace& fes,
		FaceElementTransformations& trans,
		ElementId e,
		GlobalConnectivity& map,
		GlobalPlusIsFaceNbr& plus_is_nbr);
	std::unordered_set<int> buildSharedLocalFaceSet(const ParMesh& pmesh);

}