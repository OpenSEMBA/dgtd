#include "evolution/HesthavenEvolutionMethods.h"

namespace maxwell {

	using namespace mfem;

	HesthavenFields::HesthavenFields(int size)
	{
		for (int d = X; d <= Z; d++) {
			e_[d].resize(size);
			h_[d].resize(size);
		}
	}

	FieldsInputMaps::FieldsInputMaps(const Vector& in, FiniteElementSpace& fes)
	{
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + 0 * fes.GetNDofs(), fes.GetNDofs(), 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + 1 * fes.GetNDofs(), fes.GetNDofs(), 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + 2 * fes.GetNDofs(), fes.GetNDofs(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + 3 * fes.GetNDofs(), fes.GetNDofs(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + 4 * fes.GetNDofs(), fes.GetNDofs(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + 5 * fes.GetNDofs(), fes.GetNDofs(), 1));
	}

	FieldsOutputMaps::FieldsOutputMaps(Vector& out, FiniteElementSpace& fes)
	{
		e_.push_back(Eigen::Map<Eigen::VectorXd>(out.GetData() + 0 * fes.GetNDofs(), fes.GetNDofs(), 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(out.GetData() + 1 * fes.GetNDofs(), fes.GetNDofs(), 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(out.GetData() + 2 * fes.GetNDofs(), fes.GetNDofs(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(out.GetData() + 3 * fes.GetNDofs(), fes.GetNDofs(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(out.GetData() + 4 * fes.GetNDofs(), fes.GetNDofs(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(out.GetData() + 5 * fes.GetNDofs(), fes.GetNDofs(), 1));
	}

	FieldsElementMaps::FieldsElementMaps(const Vector& in, FiniteElementSpace& fes, const ElementId& id)
	{

		Array<int> dofs;
		auto el2dofs = fes.GetElementDofs(id, dofs);

		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + id * dofs.Size() + 0 * fes.GetNDofs(), dofs.Size(), 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + id * dofs.Size() + 1 * fes.GetNDofs(), dofs.Size(), 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + id * dofs.Size() + 2 * fes.GetNDofs(), dofs.Size(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + id * dofs.Size() + 3 * fes.GetNDofs(), dofs.Size(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + id * dofs.Size() + 4 * fes.GetNDofs(), dofs.Size(), 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.GetData() + id * dofs.Size() + 5 * fes.GetNDofs(), dofs.Size(), 1));
	}

	HesthavenElementJumps::HesthavenElementJumps(HesthavenFields& in, ElementId id, Eigen::Index& elFluxSize)
	{
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.e_[X].data() + id * elFluxSize, elFluxSize, 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.e_[Y].data() + id * elFluxSize, elFluxSize, 1));
		e_.push_back(Eigen::Map<Eigen::VectorXd>(in.e_[Z].data() + id * elFluxSize, elFluxSize, 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.h_[X].data() + id * elFluxSize, elFluxSize, 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.h_[Y].data() + id * elFluxSize, elFluxSize, 1));
		h_.push_back(Eigen::Map<Eigen::VectorXd>(in.h_[Z].data() + id * elFluxSize, elFluxSize, 1));
	}

	InteriorFaceConnectivityMaps mapConnectivity(const BilinearForm* fluxMatrix)
	{
		std::vector<int> vmapM, vmapP;
		double tol{ 1e-6 };
		Array<int> cols;
		Vector vals;
		for (auto r{ 0 }; r < fluxMatrix->NumRows() / 2; r++) {
			fluxMatrix->SpMat().GetRow(r, cols, vals);
			for (auto c{ 0 }; c < cols.Size(); c++) {
				if (vals.Size() && r != cols[c] && std::abs(vals[c] - vals[cols.Find(r)]) < tol && std::abs(vals[cols.Find(r)]) > tol) {
					vmapM.emplace_back(r);
					vmapP.emplace_back(cols[c]);
				}
			}
		}
		return std::make_pair(vmapM, vmapP);
	}

	const std::vector<Source::Position> buildDoFPositions(const FiniteElementSpace& fes)
	{
		auto fec{ dynamic_cast<const L2_FECollection*>(fes.FEColl()) };
		auto vdimfes = FiniteElementSpace(fes.GetMesh(), fec, 3);
		GridFunction nodes(&vdimfes);
		fes.GetMesh()->GetNodes(nodes);
		auto dirSize{ nodes.Size() / 3 };
		std::vector<Source::Position> res;
		res.resize(dirSize);
		for (auto i{ 0 }; i < dirSize; i++) {
			res[i] = Source::Position({ nodes[i], nodes[i + dirSize], nodes[i + dirSize * 2] });
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

	void restoreOriginalAttributesAfterSubMeshing(FaceElementTransformations* faceTrans, Mesh& m, const std::map<int, Attribute>& attMap)
	{
		m.SetAttribute(faceTrans->Elem1No, attMap.at(faceTrans->Elem1No));
		if (faceTrans->Elem2No) {
			m.SetAttribute(faceTrans->Elem2No, attMap.at(faceTrans->Elem2No));
		}
	}

	void restoreOriginalAttributesAfterSubMeshing(ElementId e, Mesh& m, const std::map<int, Attribute>& attMap) 
	{
		m.SetAttribute(e, attMap.at(e));
	}

	Array<int> getFacesForElement(const Mesh& m, const ElementId e)
	{
		auto e2f = m.ElementToEdgeTable();
		e2f.Finalize();
		Array<int> res;
		e2f.GetRow(e, res);
		return res;
	}

	const int getNumFaces(const Geometry::Type& geom)
	{
		int res;
		Geometry::Dimension[geom] == 2 ? res = Geometry::NumEdges[geom] : res = Geometry::NumFaces[geom];
		return res;
	}

	const int getFaceNodeNumByGeomType(const FiniteElementSpace& fes)
	{
		int res;
		int order = fes.FEColl()->GetOrder();
		switch (fes.GetMesh()->Dimension()) {
		case 2:
			res = order + 1;
			break;
		case 3:
			fes.GetFE(0)->GetGeomType() == Geometry::Type::TETRAHEDRON ? res = ((order + 1) * (order + 2) / 2) : res = (order + 1) * (order + 1);
			break;
		default:
			throw std::runtime_error("Method only supports 2D and 3D problems.");
		}
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

	void markElementsForSubMeshing(FaceElementTransformations* faceTrans, Mesh& m)
	{
		m.SetAttribute(faceTrans->Elem1No, hesthavenMeshingTag);
		m.SetAttribute(faceTrans->Elem2No, hesthavenMeshingTag);
	}

	DynamicMatrix loadMatrixForTris(const DynamicMatrix& global, const int startRow, const int startCol)
	{
		DynamicMatrix res;
		res.resize(int(global.rows() / 2), int(global.cols() / 2));
		for (auto r{ startRow }; r < startRow + int(global.rows() / 2); r++) {
			for (auto c{ startCol }; c < startCol + int(global.cols() / 2); c++) {
				res(r - startRow, c - startCol) = global(r, c);
			}
		}
		return res;
	}

	DynamicMatrix getElemMassMatrixFromGlobal(const ElementId e, const DynamicMatrix& global, const Element::Type elType)
	{
		switch (elType) {
		case Element::Type::TRIANGLE:
			switch (e) {
			case 0:
				return loadMatrixForTris(global, 0, 0);
			case 1:
				return loadMatrixForTris(global, int(global.rows() / 2), int(global.cols() / 2));
			default:
				throw std::runtime_error("Incorrect element index for getElementMassMatrixFromGlobal");
			}
		case Element::Type::QUADRILATERAL:
			return global;
		case Element::Type::TETRAHEDRON:
			return global;
		case Element::Type::HEXAHEDRON:
			return global;
		default:
			throw std::runtime_error("Unsupported Element Type.");
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
				if (std::abs(matrix(r, c) - 0.0) >= 1e-6)
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
			bdr_marker.SetSize(fes.GetNF());
			bdr_marker = 0;
			bdr_marker[f] = 1;
			res.emplace_back(bdr_marker);
		}
		return res;
	}

	Mesh buildHesthavenRefTetrahedra()
	{
		Mesh m(3, 0, 0, 0, 3);
		m.AddVertex(0.0, 0.0, 0.0);
		m.AddVertex(2.0, 0.0, 0.0);
		m.AddVertex(0.0, 2.0, 0.0);
		m.AddVertex(0.0, 0.0, 2.0);
		m.AddTet(0, 1, 2, 3, 777);
		m.FinalizeMesh();

		return m;
	}

	std::unique_ptr<BilinearForm> assembleFaceMassBilinearForm(FiniteElementSpace& fes, Array<int>& boundaryMarker)
	{
		auto res = std::make_unique<BilinearForm>(&fes);
		res->AddBdrFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0), boundaryMarker);
		res->Assemble();
		res->Finalize();

		return res;
	}

	DynamicMatrix assembleConnectivityFaceMassMatrix(FiniteElementSpace& subFES, Array<int> boundaryMarker)
	{
		auto bf = assembleFaceMassBilinearForm(subFES, boundaryMarker);
		auto res = toEigen(*bf->SpMat().ToDenseMatrix());
		return res;
	}

	DynamicMatrix assembleEmat(FiniteElementSpace& fes, std::vector<Array<int>>& boundaryMarkers)
	{
		DynamicMatrix res;
		auto numNodesAtFace = getFaceNodeNumByGeomType(fes);
		res.resize(fes.GetNDofs(), fes.GetNF() * numNodesAtFace);
		for (auto f{ 0 }; f < fes.GetNF(); f++) {
			auto surface_matrix{ assembleConnectivityFaceMassMatrix(fes, boundaryMarkers[f]) };
			removeZeroColumns(surface_matrix);
			res.block(0, f * numNodesAtFace, fes.GetNDofs(), numNodesAtFace) = surface_matrix;
		}
		return res;
	}

	SubMesh assembleInteriorFaceSubMesh(Mesh& m, const FaceElementTransformations& trans, const FaceToGeomTag& attMap)
	{
		Array<int> volume_marker;
		volume_marker.Append(hesthavenMeshingTag);
		m.SetAttribute(trans.Elem1No, hesthavenMeshingTag);
		m.SetAttribute(trans.Elem2No, hesthavenMeshingTag);
		auto res{ SubMesh::CreateFromDomain(m, volume_marker) };
		restoreOriginalAttributesAfterSubMeshing(trans.Elem1No, m, attMap);
		restoreOriginalAttributesAfterSubMeshing(trans.Elem2No, m, attMap);
		return res;
	}

	std::unique_ptr<BilinearForm> assembleInteriorFluxMatrix(FiniteElementSpace& fes)
	{
		auto bf{ std::make_unique<BilinearForm>(&fes) };
		bf->AddInteriorFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0));
		bf->Assemble();
		bf->Finalize();
		return bf;
	}

	DynamicMatrix assembleBoundaryFluxMatrix(FiniteElementSpace& fes)
	{
		BilinearForm bf(&fes);
		bf.AddBdrFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0));
		bf.Assemble();
		bf.Finalize();
		return toEigen(*bf.SpMat().ToDenseMatrix());
	}

	void appendConnectivityMapsForInteriorFace(const FaceElementTransformations& trans, FiniteElementSpace& globalFES, FiniteElementSpace& smFES, GlobalConnectivity& map, ElementId e)
	{
		auto bf{ assembleInteriorFluxMatrix(smFES) };
		auto maps{ mapConnectivity(bf.get()) };
		Array<int> dofs1, dofs2;
		globalFES.GetElementDofs(trans.Elem1No, dofs1);
		globalFES.GetElementDofs(trans.Elem2No, dofs2);

		ConnectivityVector sortingVector;
		if (trans.Elem1No == e) {
			for (auto v{ 0 }; v < maps.first.size(); v++) {
				sortingVector.push_back(std::pair(dofs1[maps.first[v]], dofs2[maps.second[v] - dofs1.Size()]));
			}
			sort(sortingVector.begin(), sortingVector.end());
		}
		else if (trans.Elem2No == e){
			for (auto v{ 0 }; v < maps.first.size(); v++) {
				sortingVector.push_back(std::pair(dofs2[maps.second[v] - dofs1.Size()], dofs1[maps.first[v]]));
			}
			sort(sortingVector.begin(), sortingVector.end());
		}
		else {
			throw std::runtime_error("incorrect element id to element ids in FaceElementTransformation");
		}

		for (auto v{ 0 }; v < sortingVector.size(); v++) {
			map.push_back(sortingVector[v]);
		}
	}

	InteriorFaceConnectivityMaps buildConnectivityForInteriorBdrFace(const FaceElementTransformations& trans, const FiniteElementSpace& globalFES, FiniteElementSpace& smFES)
	{
		auto matrix{ assembleInteriorFluxMatrix(smFES) };
		auto maps{ mapConnectivity(matrix.get()) };
		Array<int> dofs1, dofs2;
		globalFES.GetElementDofs(trans.Elem1No, dofs1);
		globalFES.GetElementDofs(trans.Elem2No, dofs2);

		Nodes elem1, elem2;
		
		for (auto v{ 0 }; v < maps.first.size(); v++) {
			elem1.push_back(dofs1[maps.first[v]]);
			elem2.push_back(dofs2[maps.second[v] - dofs1.Size()]);
		}

		return std::make_pair(elem1, elem2);

	}

	void appendConnectivityMapsForBoundaryFace(FiniteElementSpace& globalFES, FiniteElementSpace& submeshFES, const DynamicMatrix& surfaceMatrix, GlobalConnectivity& map)
	{
		GridFunction gfParent(&globalFES);
		GridFunction gfChild(&submeshFES);
		auto mesh = dynamic_cast<SubMesh*>(submeshFES.GetMesh());
		auto transferMap = mesh->CreateTransferMap(gfChild, gfParent);
		ConstantCoefficient zero(0.0);
		double tol{ 1e-8 };
		for (auto r{ 0 }; r < surfaceMatrix.rows(); r++) {
			if (std::abs(surfaceMatrix(r, r)) > tol) {
				gfChild.ProjectCoefficient(zero);
				gfParent.ProjectCoefficient(zero);
				gfChild[r] = 1.0;
				transferMap.Transfer(gfChild, gfParent); 
				ConnectivityVector sortingVector;
				for (auto v{ 0 }; v < gfParent.Size(); v++) {
					if (gfParent[v] != 0.0) {
						sortingVector.push_back(std::make_pair(v, v));
						break;
					}
				}
				sort(sortingVector.begin(), sortingVector.end());
				for (auto v{ 0 }; v < sortingVector.size(); v++) {
					map.push_back(sortingVector[0]);
				}
			}
		}
	}

	std::vector<Array<int>> assembleBdrMarkersForBdrElements(Mesh& bdrMesh, const int numBdrElems)
	{
		std::vector<Array<int>> res(numBdrElems);
		bdrMesh.bdr_attributes.SetSize(numBdrElems);
		for (auto b{ 0 }; b < numBdrElems; b++) {
			bdrMesh.bdr_attributes[b] = b + 1;
			bdrMesh.SetBdrAttribute(b, bdrMesh.bdr_attributes[b]);
			Array<int> bdr_marker;
			bdr_marker.SetSize(numBdrElems);
			bdr_marker = 0;
			bdr_marker[b] = 1;
			res[b] = bdr_marker;
		}
		return res;
	}

	void tagBdrAttributesForSubMesh(const int edge, SubMesh& sm)
	{
		for (auto b{ 0 }; b < sm.bdr_attributes.Size(); b++) {
			sm.bdr_attributes[b] = 1;
			sm.SetBdrAttribute(b, 1);
		}
		sm.SetBdrAttribute(edge, hesthavenMeshingTag);
		sm.bdr_attributes.Append(hesthavenMeshingTag);
	}

	Nodes findNodesPerBdrFace(const BilinearForm* bdrNodeMat)
	{
		Nodes res;
		Array<int> cols;
		Vector vals;
		double tol{ 1e-6 };
		for (auto r{ 0 }; r < bdrNodeMat->NumRows(); r++) {
			bdrNodeMat->SpMat().GetRow(r, cols, vals);
			if (vals.Size() && cols.Find(r) != -1 && std::abs(vals[cols.Find(r)]) > tol) {
				res.push_back(r);
			}
		}

		return res;
	}

	std::vector<Nodes> assembleNodeVectorPerBdrFace(std::vector<Array<int>>& bdrNodeMarkers, FiniteElementSpace& bdrFES, const std::map<bool, std::vector<BdrElementId>>& isInteriorMap)
	{

		std::vector<Nodes> res(isInteriorMap.at(false).size());

		for (auto b{ 0 }; b < isInteriorMap.at(false).size(); b++) {
			auto bdrNodeFinderOperator{ assembleFaceMassBilinearForm(bdrFES, bdrNodeMarkers[isInteriorMap.at(false)[b]]) };
			res[b] = findNodesPerBdrFace(bdrNodeFinderOperator.get());
		}

		return res;
	}

	const std::map<bool, std::vector<BdrElementId>> assembleInteriorOrTrueBdrMap(const FiniteElementSpace& fes)
	{
		std::map<bool, std::vector<BdrElementId>> res;
		auto f2bdr{ fes.GetMesh()->GetFaceToBdrElMap() };
		for (auto b{ 0 }; b < fes.GetNBE(); b++) {
			fes.GetMesh()->FaceIsInterior(f2bdr.Find(b)) ==
				true ? res[true].push_back(b) : res[false].push_back(b);
		}
		return res;
	}

	Connectivities::Connectivities(Model& model, const FiniteElementSpace& fes) :
		model_(model),
		fes_(fes)
	{
		auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };
		global = assembleGlobalConnectivityMap(model.getMesh(), fec);

		const auto isInteriorMap{ assembleInteriorOrTrueBdrMap(fes_) };

		{
			auto bdrNodeMesh{ Mesh(model_.getMesh()) };
			std::vector<Array<int>> bdrNodeMarkers{ assembleBdrMarkersForBdrElements(bdrNodeMesh, model_.getConstMesh().GetNBE()) };
			auto bdrNodeFES = FiniteElementSpace(&bdrNodeMesh, fec);
			std::vector<int> atts;
			for (auto b = 0; b < model_.getConstMesh().GetNBE(); b++) {
				atts.push_back(model.getConstMesh().GetBdrAttribute(b));
			}
			const auto bdr2nodes{ assembleNodeVectorPerBdrFace(bdrNodeMarkers, bdrNodeFES, isInteriorMap) };
			initBdrConnectivityMaps(bdr2nodes, isInteriorMap);
		}

		{
			if (model_.getInteriorBoundaryToMarker().size() != 0) {
				initIntFaceConnectivityMaps(model_.getInteriorBoundaryToMarker());
			}
		}

		{
			if (model_.getTotalFieldScatteredFieldToMarker().size() != 0) {
				initIntFaceConnectivityMaps(model.getTotalFieldScatteredFieldToMarker());
			}
		}

	}

	void Connectivities::loadIntBdrConditions(const InteriorFaceConnectivityMaps& mapB, const InteriorFaceConnectivityMaps& nodePairs, const BdrCond& bdrCond, const double ori = 0.0)
	{
		if (bdrCond != BdrCond::TotalFieldIn) {
			switch (bdrCond) {
			case BdrCond::PEC:
				boundary.intPEC.mapBElem1.push_back(mapB.first);
				boundary.intPEC.mapBElem2.push_back(mapB.second);
				boundary.intPEC.vmapBElem1.push_back(nodePairs.first);
				boundary.intPEC.vmapBElem2.push_back(nodePairs.second);
				break;
			case BdrCond::PMC:
				boundary.intPMC.mapBElem1.push_back(mapB.first);
				boundary.intPMC.mapBElem2.push_back(mapB.second);
				boundary.intPMC.vmapBElem1.push_back(nodePairs.first);
				boundary.intPMC.vmapBElem2.push_back(nodePairs.second);
				break;
			case BdrCond::SMA:
				boundary.intSMA.mapBElem1.push_back(mapB.first);
				boundary.intSMA.mapBElem2.push_back(mapB.second);
				boundary.intSMA.vmapBElem1.push_back(nodePairs.first);
				boundary.intSMA.vmapBElem2.push_back(nodePairs.second);
				break;
			default:
				throw std::runtime_error("Incorrect boundary condition for interior assignment.");
			}
		}
		else if (bdrCond == BdrCond::TotalFieldIn) {
			if (ori >= 0.0) {
				boundary.TFSF.mapBSF.push_back(mapB.first);
				boundary.TFSF.mapBTF.push_back(mapB.second);
				boundary.TFSF.vmapBSF.push_back(nodePairs.first);
				boundary.TFSF.vmapBTF.push_back(nodePairs.second);
			}
			else {
				boundary.TFSF.mapBSF.push_back(mapB.second);
				boundary.TFSF.mapBTF.push_back(mapB.first);
				boundary.TFSF.vmapBSF.push_back(nodePairs.second);
				boundary.TFSF.vmapBTF.push_back(nodePairs.first);
			}
		}
	}

	InteriorFaceConnectivityMaps Connectivities::initInteriorFacesMapB(const InteriorFaceConnectivityMaps& nodePairs) const
	{
		InteriorFaceConnectivityMaps res;
		res.first.resize(nodePairs.first.size());
		res.second.resize(nodePairs.second.size());
		for (auto v{ 0 }; v < res.first.size(); v++) {
			res.first[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(nodePairs.first[v], nodePairs.second[v])));
			res.second[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(nodePairs.second[v], nodePairs.first[v])));
		}
		return res;
	}

	void Connectivities::initIntFaceConnectivityMaps(const BoundaryToMarker& markers)
	{
		Array<int> elementMarker;
		elementMarker.Append(hesthavenMeshingTag);

		auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };
		auto attMap{ mapOriginalAttributes(model_.getConstMesh()) };

		double ori{ 0.0 };
		auto intBdrNodeMesh{ Mesh(model_.getConstMesh()) };
		for (auto b{ 0 }; b < model_.getConstMesh().GetNBE(); b++) {
			for (const auto& marker : markers) {
				if (marker.second[model_.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
					auto faceTrans{ model_.getMesh().GetInternalBdrFaceTransformations(b) };
					fes_.GetMesh()->Dimension() == 2 ? ori = calculateCrossBaryVertexSign(*fes_.GetMesh(), *faceTrans, b) : ori = buildFaceOrientation(*fes_.GetMesh(), b);
					auto twoElemSubMesh{ assembleInteriorFaceSubMesh(model_.getMesh(), *faceTrans, attMap) };
					FiniteElementSpace subFES(&twoElemSubMesh, fec);
					auto nodePairs{ buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES) };
					auto mapsB{ initInteriorFacesMapB(nodePairs) };
					loadIntBdrConditions(mapsB, nodePairs, marker.first, ori);
				}
			}
		}
	}

	std::vector<NodeId> findMapBIDForFaceNodes(const Nodes& expectedNodes, const GlobalConnectivity& connectivity)
	{
		Nodes sortedNodes = expectedNodes;
		std::sort(sortedNodes.begin(), sortedNodes.end());
		std::vector<NodeId> res(sortedNodes.size());
		Nodes connNodes(sortedNodes.size());
		for (auto n{ 0 }; n < connectivity.size(); n++) {
			if (sortedNodes[0] == connectivity[n].first && n + sortedNodes.size() - 1 < connectivity.size()) {
				for (auto v{ 0 }; v < sortedNodes.size(); v++) {
					connNodes[v] = connectivity[n + v].first;
				}
				if (sortedNodes == connNodes) {
					for (auto v{ 0 }; v < sortedNodes.size(); v++) {
						res[v] = n + v;
					}
					return res;
				}
			}
		}
	}

	Nodes initMapB(const GlobalConnectivity& connectivity)
	{
		Nodes res;
		for (auto i{ 0 }; i < connectivity.size(); i++) {
			if (connectivity[i].first == connectivity[i].second) {
				res.push_back(i);
			}
		}
		return res;
	}


	Nodes initVMapB(const GlobalConnectivity& connectivity, const Nodes& mapB)
	{
		Nodes res(mapB.size());
		for (auto i{ 0 }; i < mapB.size(); i++) {
			res[i] = connectivity[mapB[i]].first;
		}
		return res;
	}

	void Connectivities::initBdrConnectivityMaps(const std::vector<Nodes>& bdr2nodes, const std::map<bool, std::vector<BdrElementId>>& isInteriorMap)
	{
		auto mapB{ initMapB(global) };
		auto vmapB{ initVMapB(global, mapB) };

		auto modelBdrMarkers = model_.getBoundaryToMarker(); //Only for true boundary types (not interior)...
		for (auto b{ 0 }; b < bdr2nodes.size(); b++) { //For each one of our bdr elements...
			for (const auto& marker : modelBdrMarkers) { //Fetch each one of our bdrMarkers in the model...
				if (marker.second[fes_.GetBdrAttribute(isInteriorMap.at(false)[b]) - 1] == 1) { //And check if that BdrCond is active for the specified bdrAtt of that bdr element, which means the bdr element is of that type...
					for (auto n{ 0 }; n < vmapB.size(); n++) { //Then sweep vmapB...
						auto mapBIds{ findMapBIDForFaceNodes(bdr2nodes[b], global) }; //To find the nodes in vmapB that are equal to the first node in our bdr element (as they are sorted when built), because it can happen...
						Nodes mapBToStore(bdr2nodes[b].size()); //... that we have multiple instances of a Hesthaven type node appearing twice on vmapB (i.e. corner node)...
						Nodes vmapBToStore(bdr2nodes[b].size()); //... so we need to save the mapB 2 vmapB pairs so the nodes are properly linked when defining boundary conditions...
						for (auto m{ 0 }; m < mapBIds.size(); m++) {
							mapBToStore[m] = mapBIds[m]; //Then make a temporary vector that we'll fill the next nodesize worth of nodes for mapB and vmapB, starting at the compared and equal node...
							vmapBToStore[m] = global[mapBIds[m]].first;
						}
						switch (marker.first) {
						case BdrCond::PEC:
							boundary.PEC.mapB.push_back(mapBToStore);
							boundary.PEC.vmapB.push_back(vmapBToStore);
							break;
						case BdrCond::PMC:
							boundary.PMC.mapB.push_back(mapBToStore);
							boundary.PMC.vmapB.push_back(vmapBToStore);
							break;
						case BdrCond::SMA:
							boundary.SMA.mapB.push_back(mapBToStore);
							boundary.SMA.vmapB.push_back(vmapBToStore);
							break;
						}
						break; // And as we only support one bdrcond per face/edge/bdrwhatever, we're done for this bdr element...
					}     // As this method is completely non-dependent on dimensions or geometries, could be assumed is as generic as it could get.
				}
			}
		}
	}

	GlobalConnectivity assembleGlobalConnectivityMap(Mesh& m, const L2_FECollection* fec)
	{
		GlobalConnectivity res;
		auto mesh{ Mesh(m) };
		FiniteElementSpace globalFES(&mesh, fec);

		std::map<FaceId, bool> global_face_is_interior;
		int numFaces;
		m.Dimension() == 2 ? numFaces = mesh.GetNEdges() : numFaces = mesh.GetNFaces();
		for (auto f{ 0 }; f < numFaces; f++) {
			global_face_is_interior[f] = mesh.FaceIsInterior(f);
		}

		Table globalElementToFace;
		m.Dimension() == 2 ? globalElementToFace = mesh.ElementToEdgeTable() : globalElementToFace = mesh.ElementToFaceTable();

		Array<int> volumeMarker;
		volumeMarker.Append(hesthavenMeshingTag);

		Array<int> boundaryMarker(hesthavenMeshingTag);
		boundaryMarker = 0;
		boundaryMarker[hesthavenMeshingTag - 1] = 1;

		auto attMap{ mapOriginalAttributes(mesh) };


		for (auto e{ 0 }; e < mesh.GetNE(); e++) {

			Array<int> localFaceIndexToGlobalFaceIndex;
			globalElementToFace.GetRow(e, localFaceIndexToGlobalFaceIndex);

			int numLocalFaces;
			mesh.Dimension() == 2 ? numLocalFaces = mesh.GetElement(e)->GetNEdges() : numLocalFaces = mesh.GetElement(e)->GetNFaces();
			for (auto localFace{ 0 }; localFace < numLocalFaces; localFace++) {

				if (!global_face_is_interior[localFaceIndexToGlobalFaceIndex[localFace]]) {

					mesh.SetAttribute(e, hesthavenMeshingTag);
					auto sm = SubMesh::CreateFromDomain(mesh, volumeMarker);
					restoreOriginalAttributesAfterSubMeshing(e, mesh, attMap);
					tagBdrAttributesForSubMesh(localFace, sm);
					FiniteElementSpace smFES(&sm, fec);
					appendConnectivityMapsForBoundaryFace(
						globalFES,
						smFES,
						assembleConnectivityFaceMassMatrix(smFES, boundaryMarker),
						res);

				}
				else {

					FaceElementTransformations* faceTrans = mesh.GetFaceElementTransformations(localFaceIndexToGlobalFaceIndex[localFace]);
					auto sm = assembleInteriorFaceSubMesh(mesh, *faceTrans, attMap);
					restoreOriginalAttributesAfterSubMeshing(faceTrans, mesh, attMap);
					FiniteElementSpace smFES(&sm, fec);
					appendConnectivityMapsForInteriorFace(
						*faceTrans,
						globalFES,
						smFES,
						res,
						e);

				}
			}
		}
		return res;
	}
}