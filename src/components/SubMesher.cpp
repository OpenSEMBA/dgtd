#include "SubMesher.h"

#include <mfem.hpp>



namespace maxwell {

using namespace mfem;

TotalFieldScatteredField_SubMesher::TotalFieldScatteredField_SubMesher(const Mesh& m)
{
	Mesh parent(m);
	setTFSFAttributesForSubMeshing(parent);

	Array<int> tf_att(1); tf_att[0] = 1000;
	Array<int> sf_att(1); sf_att[0] = 2000;
	auto tf_sm{ SubMesh::CreateFromDomain(parent, tf_att) };
	auto sf_sm{ SubMesh::CreateFromDomain(parent, sf_att) };

	auto tf_m2sm_map{ SubMeshUtils::BuildFaceMap(parent, tf_sm, tf_sm.GetParentElementIDMap()) };
	auto sf_m2sm_map{ SubMeshUtils::BuildFaceMap(parent, sf_sm, sf_sm.GetParentElementIDMap()) };
	auto f2bdr_map{ parent.GetFaceToBdrElMap() };
	for (int i = 0; i < parent.GetNBE(); i++) {
		if (parent.GetBdrAttribute(i) == 301) {
			tf_sm.SetBdrAttribute(tf_m2sm_map.Find(f2bdr_map.Find(i)), 301);
			sf_sm.SetBdrAttribute(sf_m2sm_map.Find(f2bdr_map.Find(i)), 301);
		}
	}
	tf_sm.FinalizeMesh();
	sf_sm.FinalizeMesh();
};


void TotalFieldScatteredField_SubMesher::setAttributeForTagging(Mesh& m, const FaceElementTransformations* trans, bool el1_is_tf = true)
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

void TotalFieldScatteredField_SubMesher::storeElementToFaceInformation(const FaceElementTransformations* trans, const FaceId f, bool el1_is_tf = true)
{
	switch (el1_is_tf) {
	case true:
		elem_to_face_tf_.push_back(std::make_pair(trans->Elem1No, f));
		elem_to_face_sf_.push_back(std::make_pair(trans->Elem2No, (f + 2) % 4));
		break;
	case false:
		elem_to_face_sf_.push_back(std::make_pair(trans->Elem1No, (f + 2) % 4));
		elem_to_face_tf_.push_back(std::make_pair(trans->Elem2No, f));
		break;
	}
}

void TotalFieldScatteredField_SubMesher::prepareSubMeshInfo(Mesh& m, const FaceElementTransformations* trans, const FaceId f, bool el1_is_tf = true)
{
	setAttributeForTagging(m, trans, el1_is_tf);
	storeElementToFaceInformation(trans, f, el1_is_tf);
}

void TotalFieldScatteredField_SubMesher::setTFSFAttributesForSubMeshing(Mesh& m)
{

	for (int be = 0; be < m.GetNBE(); be++)
	{
		if (m.GetBdrAttribute(be) == 301)
		{
			Array<int> be_vert(2);
			Array<int> el1_vert(4);
			Array<int> el2_vert(4);

			auto be_trans{ m.GetInternalBdrFaceTransformations(be) };
			auto el1{ m.GetElement(be_trans->Elem1No) };
			auto el2{ m.GetElement(be_trans->Elem2No) };

			auto ver1{ el1->GetVertices() };
			auto ver2{ el2->GetVertices() };
			m.GetBdrElementVertices(be, be_vert);

			for (int v = 0; v < el1_vert.Size(); v++) {
				el1_vert[v] = ver1[v];
				el2_vert[v] = ver2[v];
			}

			std::vector<El2Att> el_to_att_pre(2);

			//be_vert is counterclockwise, that is our convention to designate which element will be TF. The other element will be SF.

			for (int v = 0; v < el1_vert.Size(); v++) {

				if (el1->GetAttribute() == 1000 && el2->GetAttribute() == 2000 || el1->GetAttribute() == 2000 && el2->GetAttribute() == 1000) {
					continue;
				}

				if (v + 1 < el1_vert.Size() && el1_vert[v] == be_vert[0] && el1_vert[v + 1] == be_vert[1]) {
					prepareSubMeshInfo(m, be_trans, v);
				}
				else if (v + 1 < el2_vert.Size() && el2_vert[v] == be_vert[0] && el2_vert[v + 1] == be_vert[1]) {
					prepareSubMeshInfo(m, be_trans, v, false);
				}

				//It can happen the vertex to check will not be adjacent in the Array but will be in the last and first position, as if closing the loop, thus we require an extra check.
				if (v == el1_vert.Size() - 1 && el1->GetAttribute() != 1000 && el2->GetAttribute() != 2000 || v == el1_vert.Size() - 1 && el1->GetAttribute() != 2000 && el2->GetAttribute() != 1000) {
					if (el1_vert[el1_vert.Size() - 1] == be_vert[0] && el1_vert[0] == be_vert[1]) {
						prepareSubMeshInfo(m, be_trans, v);
					}
					else if (el2_vert[el2_vert.Size() - 1] == be_vert[0] && el2_vert[0] == be_vert[1]) {
						prepareSubMeshInfo(m, be_trans, v, false);
					}
				}
			}
		}
	}
}
};