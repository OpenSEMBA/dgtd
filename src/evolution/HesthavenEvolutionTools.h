#pragma once

#include <unordered_set>

#include "mfemExtension/BilinearIntegrators.h"
#include "math/EigenMfemTools.h"
#include "components/Types.h"
#include "components/Model.h"

namespace maxwell {

	using namespace mfem;

	using NodeId = int;
	using Volume = double;
	using Nodes = std::vector<NodeId>;
	using DynamicMatrix = Eigen::MatrixXd;
	using InteriorFaceConnectivityMaps = std::pair<Nodes, Nodes>;
	using NodePair = std::pair<NodeId, NodeId>;
	using ConnectivityVector = std::vector<NodePair>;
	using GlobalConnectivity = ConnectivityVector;
	using BdrCondToNodes = std::pair<BdrCond, Nodes>;
	using Emat = std::vector<const DynamicMatrix*>;
	using Directional = std::array<const DynamicMatrix*, 3>;
	using Normals = std::array<Eigen::VectorXd, 3>;
	using TFSFSigns = std::pair<double, double>;
	using TFSFConnectivity = std::vector<std::pair<NodePair, TFSFSigns>>;

	struct MatrixCompareLessThan {

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
		Volume vol;
		Directional der;
		Emat emat;
		const DynamicMatrix* invmass;
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
		std::vector<Eigen::Map<Eigen::VectorXd>> e_, h_;
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
	void applyBoundaryConditionsToNodes(const BoundaryMaps&, const FieldsInputMaps& in, HesthavenFields& out);
	void applyInteriorBoundaryConditionsToNodes(const InteriorBoundaryMaps&, const FieldsInputMaps& in, HesthavenFields& out);

	const int getNumFaces(const Geometry::Type&); 
	const int getFaceNodeNumByGeomType(const FiniteElementSpace& fes);

	std::map<int, Attribute> mapOriginalAttributes(const Mesh& m);

	DynamicMatrix loadMatrixWithValues(const DynamicMatrix& global, const int startRow, const int startCol);
	DynamicMatrix getElementMassMatrixFromGlobal(const int el, const DynamicMatrix& global);
	DynamicMatrix getFaceMassMatrixFromGlobal(const DynamicMatrix& global);
	DynamicMatrix assembleConnectivityFaceMassMatrix(FiniteElementSpace&, Array<int> boundaryMarker);
	DynamicMatrix assembleEmat(FiniteElementSpace&, std::vector<Array<int>>& boundaryMarkers);
	std::unique_ptr<BilinearForm> assembleInteriorFluxMatrix(FiniteElementSpace&);
	DynamicMatrix assembleBoundaryFluxMatrix(FiniteElementSpace&);
	
	SubMesh assembleInteriorFaceSubMesh(Mesh& m, const FaceElementTransformations& trans, const FaceToGeomTag& attMap);
	
	std::pair<Nodes, Nodes> buildConnectivityForInteriorBdrFace(const FaceElementTransformations&, FiniteElementSpace& globalFES, FiniteElementSpace& smFES);
	std::vector<Array<int>> assembleBoundaryMarkers(FiniteElementSpace&);
	std::unique_ptr<BilinearForm> assembleFaceMassBilinearForm(FiniteElementSpace&, Array<int>& boundaryMarker);
	InteriorFaceConnectivityMaps mapConnectivity(const BilinearForm* fluxMatrix);
	Array<int> getFacesForElement(const Mesh&, const ElementId);
	FaceElementTransformations* getInteriorFaceTransformation(Mesh&, const Array<int>& faces);
	GlobalConnectivity assembleGlobalConnectivityMap(Mesh&, const L2_FECollection*);

}