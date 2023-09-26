#include "SubMesher.h"

#include <mfem.hpp>

namespace maxwell {

using namespace mfem;

TotalFieldScatteredFieldSubMesher::TotalFieldScatteredFieldSubMesher(const Mesh& m)
{
	Mesh parent(m);
 	setTFSFAttributesForSubMeshing(parent);

	Array<int> tf_att(1); tf_att[0] = 1000;
	Array<int> sf_att(1); sf_att[0] = 2000;
	auto tf_sm{ SubMesh::CreateFromDomain(parent, tf_att) };
	auto sf_sm{ SubMesh::CreateFromDomain(parent, sf_att) };

	auto f2bdr_map{ parent.GetFaceToBdrElMap() };
	
	setBoundaryAttributesInChild(parent,tf_sm);
	setBoundaryAttributesInChild(parent,sf_sm);
	
	restoreElementAttributes(tf_sm);
	restoreElementAttributes(sf_sm);

	tf_sm.FinalizeMesh();
	sf_sm.FinalizeMesh();

	tf_mesh_ = std::move(tf_sm);
	sf_mesh_ = std::move(sf_sm);

};

void TotalFieldScatteredFieldSubMesher::restoreElementAttributes(Mesh& m) //Temporary method that has to be reworked when materials are implemented.
{
	for (int e = 0; e < m.GetNE(); e++)
	{
		m.SetAttribute(e, 1);
	}
}

void TotalFieldScatteredFieldSubMesher::setBoundaryAttributesInChild(const Mesh& parent, SubMesh& child)
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
	switch (el1_is_tf) {
	case true:
		m.GetElement(trans->Elem1No)->SetAttribute(1000);
		m.GetElement(trans->Elem2No)->SetAttribute(2000);
		break;
	case false:
		m.GetElement(trans->Elem1No)->SetAttribute(2000);
		m.GetElement(trans->Elem2No)->SetAttribute(1000);
		break;
	}
}

void TotalFieldScatteredFieldSubMesher::storeElementToFaceInformation(const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf)
{
	switch (el1_is_tf) {
	case true:
		elem_to_face_tf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		elem_to_face_sf_.push_back(std::make_pair(trans->Elem2No, facesId.second));
		break;
	case false:
		elem_to_face_tf_.push_back(std::make_pair(trans->Elem2No, facesId.second));
		elem_to_face_sf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
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
		if (v == verts.Size() - 1 && el->GetAttribute() != 1000 || v == verts.Size() - 1 && el->GetAttribute() != 2000) {
			if (verts[verts.Size() - 1] == be_verts[0] && verts[0] == be_verts[1]) {
				return std::make_pair(v, true);
			}
		}

		//Check for clockwise.
		if (v + 1 < verts.Size() && verts[v] == be_verts[1] && verts[v + 1] == be_verts[0]) {
			return std::make_pair(v, false);
		}

		//It can happen the vertices to check will not be adjacent in the Array but will be in the last and first position, as if closing the loop, thus we require an extra check.
		if (v == verts.Size() - 1 && el->GetAttribute() != 1000 || v == verts.Size() - 1 && el->GetAttribute() != 2000) {
			if (verts[verts.Size() - 1] == be_verts[1] && verts[0] == be_verts[0]) {
				return std::make_pair(v, false);
			}
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setTFSFAttributesForSubMeshing(Mesh& m)
{
	for (int be = 0; be < m.GetNBE(); be++)
	{
		if (m.GetBdrAttribute(be) == 301)
		{

			auto be_trans{ m.GetInternalBdrFaceTransformations(be) };
			auto el1{ m.GetElement(be_trans->Elem1No) };
			auto el2{ m.GetElement(be_trans->Elem2No) };

			Array<int> el1_vert(el1->GetNVertices());
			Array<int> el2_vert(el2->GetNVertices());
			Array<int> be_vert(2);

			auto ver1{ el1->GetVertices() };
			auto ver2{ el2->GetVertices() };
			m.GetBdrElementVertices(be, be_vert);

			for (int v = 0; v < el1_vert.Size(); v++) {
				el1_vert[v] = ver1[v];
			}
			for (int v2 = 0; v2 < el2_vert.Size(); v2++) {
				el2_vert[v2] = ver2[v2];
			}

			//be_vert is counterclockwise, that is our convention to designate which element will be TF. The other element will be SF.

			auto set_v1 = getFaceAndDirOnVertexIteration(el1, el1_vert, be_vert);
			auto set_v2 = getFaceAndDirOnVertexIteration(el2, el2_vert, be_vert);
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

};