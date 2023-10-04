#include "SubMesher.h"

namespace maxwell {

using namespace mfem;
static const int NotFound{ -1 };

TotalFieldScatteredFieldSubMesher::TotalFieldScatteredFieldSubMesher(const Mesh& m)
{
	Mesh parent_for_global(m);
	Mesh parent_for_individual(m);

	setGlobalTFSFAttributesForSubMeshing(parent_for_global);
	switch (parent_for_individual.Dimension()) {
	case 1:
		setIndividualTFSFAttributesForSubMeshing1D(parent_for_individual);
		break;
	case 2:
		setIndividualTFSFAttributesForSubMeshing2D(parent_for_individual);
		break;
	}

	Array<int> global_att(1); global_att[0] = SubMeshingMarkers::Global_SubMesh;
	auto global_sm{ SubMesh::CreateFromDomain(parent_for_global, global_att) };
	restoreElementAttributes(global_sm);
	global_sm.FinalizeMesh();
	global_submesh_ = std::make_unique<SubMesh>(global_sm);

	if (!elem_to_face_tf_.empty()) {
		tf_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(parent_for_individual, true));
	}

	if (!elem_to_face_sf_.empty()) {
		sf_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(parent_for_individual, false));
	}

};

SubMesh TotalFieldScatteredFieldSubMesher::createSubMeshFromParent(const Mesh& parent, bool isTF)
{
	Array<int> marker(1); 
	if (isTF) {
		marker[0] = SubMeshingMarkers::TotalField;
	}
	else {
		marker[0] = SubMeshingMarkers::ScatteredField;
	}
	auto res{ SubMesh::CreateFromDomain(parent, marker) };
	switch (parent.Dimension()) {
	case 1:
		setBoundaryAttributesInChild1D(parent, res);
		break;
	case 2:
		setBoundaryAttributesInChild2D(parent, res);
		break;
	}
	restoreElementAttributes(res);
	res.FinalizeMesh();
	return res;
}

void TotalFieldScatteredFieldSubMesher::restoreElementAttributes(Mesh& m) //Temporary method that has to be reworked when materials are implemented.
{
	for (int e = 0; e < m.GetNE(); e++)
	{
		m.SetAttribute(e, 1);
	}
}

void TotalFieldScatteredFieldSubMesher::setBoundaryAttributesInChild1D(const Mesh& parent, SubMesh& child)
{
	for (int e = 0; e < child.GetNE(); e++) {
		Array<int> verts(2);
		child.GetElementVertices(e, verts);
		child.AddBdrPoint(verts[0]);
		child.AddBdrPoint(verts[1]);
	}
	auto parent_f2bdr_map{ parent.GetFaceToBdrElMap() };
	auto child_f2bdr_map{ child.GetFaceToBdrElMap() };
	auto map{ SubMeshUtils::BuildFaceMap(parent, child, child.GetParentElementIDMap()) };
	for (int i = 0; i < parent.GetNBE(); i++) {
		if (parent.GetBdrAttribute(i) == 301) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(i))], 301);
		}
		else if (parent.GetBdrAttribute(i) == 302) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(i))], 302);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setBoundaryAttributesInChild2D(const Mesh& parent, SubMesh& child)
{
	auto parent_f2bdr_map{ parent.GetFaceToBdrElMap() };
	auto child_f2bdr_map{ child.GetFaceToBdrElMap() };
	auto map{ SubMeshUtils::BuildFaceMap(parent, child, child.GetParentElementIDMap()) };
	for (int i = 0; i < parent.GetNBE(); i++) {
		if (parent.GetBdrAttribute(i) == 301) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(i))], 301);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setAttributeForTagging(Mesh& m, const FaceElementTransformations* trans, bool el1_is_tf)
{
	if (trans->Elem2No != NotFound) {
		if (el1_is_tf) {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::TotalField);
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::ScatteredField);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::ScatteredField);
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::TotalField);
		}
	}
	else {
		if (el1_is_tf) {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::TotalField);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::ScatteredField);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::storeElementToFaceInformation(const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf)
{
	if (facesId.second != NotFound) {
		if (el1_is_tf) {
			elem_to_face_tf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
			elem_to_face_sf_.push_back(std::make_pair(trans->Elem2No, facesId.second));
		}
		else {
			elem_to_face_tf_.push_back(std::make_pair(trans->Elem2No, facesId.second));
			elem_to_face_sf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		}
	}
	else {
		if (el1_is_tf) {
			elem_to_face_tf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		}
		else {
			elem_to_face_sf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		}
	}
}

void TotalFieldScatteredFieldSubMesher::prepareSubMeshInfo(Mesh& m, const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf)
{
	setAttributeForTagging(m, trans, el1_is_tf);
	storeElementToFaceInformation(trans, facesId, el1_is_tf);
}

Face2Dir TotalFieldScatteredFieldSubMesher::getFaceAndDirOnVertexIteration(const Element* el, const Array<int>& verts, const Array<int>& be_verts)
{
	for (int v = 0; v < verts.Size(); v++) {

		//Check for counterclockwise.
		if (v + 1 < verts.Size() && verts[v] == be_verts[0] && verts[v + 1] == be_verts[1]) {
			return std::make_pair(v, true);
		}

		//It can happen the vertices to check will not be adjacent in the Array but will be in the last and first position, as if closing the loop, thus we require an extra check.
		if (v == verts.Size() - 1 && el->GetAttribute() != SubMeshingMarkers::TotalField || v == verts.Size() - 1 && el->GetAttribute() != SubMeshingMarkers::ScatteredField) {
			if (verts[verts.Size() - 1] == be_verts[0] && verts[0] == be_verts[1]) {
				return std::make_pair(v, true);
			}
		}

		//Check for clockwise.
		if (v + 1 < verts.Size() && verts[v] == be_verts[1] && verts[v + 1] == be_verts[0]) {
			return std::make_pair(v, false);
		}

		//It can happen the vertices to check will not be adjacent in the Array but will be in the last and first position, as if closing the loop, thus we require an extra check.
		if (v == verts.Size() - 1 && el->GetAttribute() != SubMeshingMarkers::TotalField || v == verts.Size() - 1 && el->GetAttribute() != SubMeshingMarkers::ScatteredField) {
			if (verts[verts.Size() - 1] == be_verts[1] && verts[0] == be_verts[0]) {
				return std::make_pair(v, false);
			}
		}
	}
}

FaceElementTransformations* TotalFieldScatteredFieldSubMesher::getFaceElementTransformation(Mesh&m, int be) 
{
	switch (m.FaceIsInterior(m.GetBdrFace(be))) {
	case true:
		return m.GetInternalBdrFaceTransformations(be);
	default:
		return m.GetBdrFaceTransformations(be);
	}
}

void TotalFieldScatteredFieldSubMesher::setGlobalTFSFAttributesForSubMeshing(Mesh& m)
{

	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn) || 
			m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {
			auto be_trans{ getFaceElementTransformation(m, be)};
			if (be_trans->Elem2No != NotFound) {
				m.GetElement(be_trans->Elem1No)->SetAttribute(SubMeshingMarkers::Global_SubMesh);
				m.GetElement(be_trans->Elem2No)->SetAttribute(SubMeshingMarkers::Global_SubMesh);
				elems_for_global_submesh_.push_back(be_trans->Elem1No);
				elems_for_global_submesh_.push_back(be_trans->Elem2No);
			}
			else {
				m.GetElement(be_trans->Elem1No)->SetAttribute(SubMeshingMarkers::Global_SubMesh);
				elems_for_global_submesh_.push_back(be_trans->Elem1No);
			}
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing1D(Mesh& m)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn) || 
			m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {
			auto be_trans{ getFaceElementTransformation(m, be) };
			auto el1{ m.GetElement(be_trans->Elem1No) };

			auto ver1{ el1->GetVertices() };
			Array<int> be_vert(1);
			m.GetBdrElementVertices(be, be_vert);
			if (m.GetElement(be_trans->Elem2No) != NULL) {

				auto el2{ m.GetElement(be_trans->Elem2No) };
				auto ver2{ el2->GetVertices() };

				if (ver1[1] == be_vert[0] && ver2[0] == be_vert[0] && m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn)) {
					el1->SetAttribute(SubMeshingMarkers::ScatteredField);
					el2->SetAttribute(SubMeshingMarkers::TotalField);
					elem_to_face_tf_.push_back(std::make_pair(be_trans->Elem2No, ver2[0]));
					elem_to_face_sf_.push_back(std::make_pair(be_trans->Elem1No, ver1[1]));
				}
				if (ver1[1] == be_vert[0] && ver2[0] == be_vert[0] && m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {
					el1->SetAttribute(SubMeshingMarkers::TotalField);
					el2->SetAttribute(SubMeshingMarkers::ScatteredField);
					elem_to_face_tf_.push_back(std::make_pair(be_trans->Elem1No, ver1[1]));
					elem_to_face_sf_.push_back(std::make_pair(be_trans->Elem2No, ver2[0]));
				}
			}
			else {

				if (ver1[1] == be_vert[0] && m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn)) {
					el1->SetAttribute(SubMeshingMarkers::ScatteredField);
					elem_to_face_sf_.push_back(std::make_pair(be_trans->Elem1No, ver1[1]));
				}
				if (ver1[1] == be_vert[0] && m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {
					el1->SetAttribute(SubMeshingMarkers::TotalField);
					elem_to_face_tf_.push_back(std::make_pair(be_trans->Elem1No, ver1[1]));
				}
			}
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing2D(Mesh& m)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn)) {

			auto be_trans{ getFaceElementTransformation(m, be) };
			auto el1{ m.GetElement(be_trans->Elem1No) };

			Array<int> el1_vert(el1->GetNVertices());
			{
				auto ver1{ el1->GetVertices() };
				for (int v = 0; v < el1_vert.Size(); v++) {
					el1_vert[v] = ver1[v];
				}
			}
			Array<int> be_vert(2);
			m.GetBdrElementVertices(be, be_vert);
			
			auto set_v1 = getFaceAndDirOnVertexIteration(el1, el1_vert, be_vert);
			Face2Dir set_v2;

			if (be_trans->Elem2No != NotFound) {

				auto el2{ m.GetElement(be_trans->Elem2No) };

				Array<int> el2_vert(el2->GetNVertices());
				{
					auto ver2{ el2->GetVertices() };
					for (int v2 = 0; v2 < el2_vert.Size(); v2++) {
						el2_vert[v2] = ver2[v2];
					}
				}

				set_v2 = getFaceAndDirOnVertexIteration(el2, el2_vert, be_vert);
			}
			else {
				set_v2 = std::make_pair(NotFound, false);
			}
			//be_vert is counterclockwise, that is our convention to designate which element will be TF. The other element will be SF.
			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_v1.first, set_v2.first);
			std::pair<IsCCW, IsCCW> dirInfo = std::make_pair(set_v1.second, set_v2.second);
			prepareSubMeshInfo(m, be_trans, facesInfo, set_v1.second);
		}
	}
}


MaxwellTransferMap::MaxwellTransferMap(const GridFunction& src,
	const GridFunction& dst)
{
	const FiniteElementSpace* parentfes = nullptr, * subfes = nullptr;
	SubMesh* src_sm = static_cast<SubMesh*>(src.FESpace()->GetMesh());
	subfes = src.FESpace();
	parentfes = dst.FESpace();
	SubMeshUtils::BuildVdofToVdofMap(*subfes,
		*parentfes,
		src_sm->GetFrom(),
		src_sm->GetParentElementIDMap(),
		sub_to_parent_map_);
}

void MaxwellTransferMap::TransferAdd(const GridFunction& src, GridFunction& dst) const
{
	for (int i = 0; i < sub_to_parent_map_.Size(); i++)
	{
		dst(sub_to_parent_map_[i]) += src(i);
	}
}

void MaxwellTransferMap::TransferSub(const GridFunction& src, GridFunction& dst) const
{
	for (int i = 0; i < sub_to_parent_map_.Size(); i++)
	{
		dst(sub_to_parent_map_[i]) -= src(i);
	}
}

};