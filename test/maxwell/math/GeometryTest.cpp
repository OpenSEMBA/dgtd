#include <gtest/gtest.h>
#include <math.h>

#include "TestUtils.h"
#include "math/Geometry.h"

using namespace maxwell;
using namespace mfem;

using FaceId = int;
using ElementId = int;
using Attribute = int;
using ElNo2Att = std::pair<ElementId, Attribute>;

class GeometryTest : public ::testing::Test {
protected:

	void backupElementAttributesPreTag(const Mesh& m, const FaceElementTransformations* trans, bool el1_is_tf = true)
	{
		switch (el1_is_tf) {
		case true:
			if (m.GetElement(trans->Elem1No)->GetAttribute() != 1000) {
				elem_to_att_tf.push_back(std::make_pair(trans->Elem1No, m.GetElement(trans->Elem1No)->GetAttribute()));
			}
			if (m.GetElement(trans->Elem2No)->GetAttribute() != 2000) {
				elem_to_att_sf.push_back(std::make_pair(trans->Elem2No, m.GetElement(trans->Elem2No)->GetAttribute()));
			}
			break;
		case false:
			if (m.GetElement(trans->Elem1No)->GetAttribute() != 2000) {
				elem_to_att_sf.push_back(std::make_pair(trans->Elem1No, m.GetElement(trans->Elem1No)->GetAttribute()));
			}
			if (m.GetElement(trans->Elem2No)->GetAttribute() != 1000) {
				elem_to_att_tf.push_back(std::make_pair(trans->Elem2No, m.GetElement(trans->Elem2No)->GetAttribute()));
			}
			break;
		}
	}

	void setAttributeForTagging(Mesh& m, const FaceElementTransformations* trans, bool el1_is_tf = true)
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

	void storeElementToFaceInformation(const FaceElementTransformations* trans, const FaceId f, bool el1_is_tf = true)
	{
		switch (el1_is_tf) {
		case true:
			elem_to_face_tf.push_back(std::make_pair(trans->Elem1No, f));
			elem_to_face_sf.push_back(std::make_pair(trans->Elem2No, (f + 2) % 4));
			break;
		case false:
			elem_to_face_sf.push_back(std::make_pair(trans->Elem1No, (f + 2) % 4));
			elem_to_face_tf.push_back(std::make_pair(trans->Elem2No, f));
			break;
		}
	}

	void prepareSubMeshInfo(Mesh& m, const FaceElementTransformations* trans, const FaceId f, bool el1_is_tf = true)
	{
		backupElementAttributesPreTag(m, trans, el1_is_tf);
		setAttributeForTagging(m, trans, el1_is_tf);
		storeElementToFaceInformation(trans, f, el1_is_tf);
	}

	void setTFSFAttributesForSubMeshing(Mesh& m) 
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

				std::vector<ElNo2Att> el_to_att_pre(2);

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

	void restoreOriginalAttributes(Mesh& m)
	{
		for (int e = 0; e < elem_to_att_tf.size(); e++)
		{
			m.SetAttribute(elem_to_att_tf[e].first, elem_to_att_tf[e].second);
		}
		for (int e = 0; e < elem_to_att_sf.size(); e++)
		{
			m.SetAttribute(elem_to_att_sf[e].first, elem_to_att_sf[e].second);
		}
	}

	void checkIfElementsHaveAttribute(const Mesh& m, const std::vector<int>& elems, const int att)
	{
		for (int e = 0; e < elems.size(); e++) {
			EXPECT_TRUE(m.GetAttribute(elems[e]) == att);
		}
	}

	std::vector<std::pair<ElementId, FaceId>> elem_to_face_tf;
	std::vector<std::pair<ElementId, FaceId>> elem_to_face_sf;
	std::vector<std::pair<ElementId, Attribute>> elem_to_att_tf;
	std::vector<std::pair<ElementId, Attribute>> elem_to_att_sf;
};

TEST_F(GeometryTest, pointsOrientation)
{
	int v1{ 1 };
	Point p1(&v1);

	ASSERT_ANY_THROW(elementsHaveSameOrientation(&p1, &p1));
}

TEST_F(GeometryTest, segmentsOrientation) 
{
	Segment s1(1, 4, 20), s2(4, 1, 30);

	EXPECT_FALSE(elementsHaveSameOrientation(&s1, &s2));
	EXPECT_TRUE(elementsHaveSameOrientation(&s1, &s1));
}

TEST_F(GeometryTest, trianglesOrientation)
{
	Triangle t1(1, 2, 3);
	Triangle t2(3, 1, 2);
	Triangle t3(3, 2, 1);

	EXPECT_TRUE(elementsHaveSameOrientation(&t1, &t1));
	EXPECT_TRUE(elementsHaveSameOrientation(&t1, &t2));
	EXPECT_FALSE(elementsHaveSameOrientation(&t1, &t3));
}

TEST_F(GeometryTest, internal_boundary_orientation_2D)
{
	//      f2         f6
	//   V3 ----- V4 ------ V5
	//   |        |         |
	// f3|  El0   |f1  El1  |f5
	//   |        |         |
	//   V0 ----- V1 ------ V2
	//      f0         f4 
	auto m{ Mesh::MakeCartesian2D(2, 1, Element::Type::QUADRILATERAL) };
	Segment seg1(1, 4, 20);
	Segment seg2(4, 1, 30);

	int bdr1{ m.AddBdrElement(new Segment(seg1)) };
	int bdr2{ m.AddBdrElement(new Segment(seg2)) };
	m.FinalizeTopology();
	m.Finalize();

	EXPECT_TRUE(elementsHaveSameOrientation(&seg1, m.GetFace(m.GetBdrFace(bdr1))));
	EXPECT_FALSE(elementsHaveSameOrientation(&seg2, m.GetFace(m.GetBdrFace(bdr2))));

}

TEST_F(GeometryTest, orientation_from_gmsh_mesh)
{
	auto m{ mfem::Mesh::LoadFromFile((gmshMeshesFolder() + "twosquares.msh").c_str(), 1, 0, true) };

	EXPECT_TRUE(elementsHaveSameOrientation(m.GetFace( m.GetBdrFace(0)), m.GetFace(m.GetBdrFace(1) )));

	EXPECT_FALSE(elementsHaveSameOrientation(m.GetBdrElement(0), m.GetBdrElement(1)));
}

TEST_F(GeometryTest, marking_element_att_through_boundary_2D)
{
	
	{
		auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square3x3marked.mesh").c_str(), 1, 0, true) };

		setTFSFAttributesForSubMeshing(m);

		checkIfElementsHaveAttribute(m, std::vector<int>{ {(1, 3, 5, 7)} }, 2000);
		checkIfElementsHaveAttribute(m, std::vector<int>{ {(4)} }, 1000);
	}

	{
		auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square5x5marked.mesh").c_str(), 1, 0, true) };

		setTFSFAttributesForSubMeshing(m);

		checkIfElementsHaveAttribute(m, std::vector<int>{ {(1, 2, 3, 5, 9, 10, 14, 15, 19, 21, 22, 23)} }, 2000);
		checkIfElementsHaveAttribute(m, std::vector<int>{ {(6, 7, 8, 11, 13, 16, 17, 18)} }, 1000);
	}

}

TEST_F(GeometryTest, marking_element_to_face_pairs_for_submeshing_2D)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square3x3marked.mesh").c_str(), 1, 0, true) };

	setTFSFAttributesForSubMeshing(m);

	std::vector<std::pair<ElementId, FaceId>>elem_to_face_tf_check{ {{4,0},{4,1},{4,2},{4,3}} };
	std::vector<std::pair<ElementId, FaceId>>elem_to_face_sf_check{ {{3,2},{7,3},{5,0},{1,1}} };

	for (int p = 0; p < elem_to_face_sf.size(); p++) {
		EXPECT_EQ(elem_to_face_tf_check, elem_to_face_tf);
		EXPECT_EQ(elem_to_face_sf_check, elem_to_face_sf);
	}

}

TEST_F(GeometryTest, transfer_parent_bdr_att_to_child)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square3x3marked.mesh").c_str(), 1, 0, true) };

	setTFSFAttributesForSubMeshing(m);

	Array<int> tf_att(1); tf_att[0] = 1000;
	Array<int> sf_att(1); sf_att[0] = 2000;
	auto tf_sm{ SubMesh::CreateFromDomain(m, tf_att) };
	auto sf_sm{ SubMesh::CreateFromDomain(m, sf_att) };

	restoreOriginalAttributes(m);

	auto tf_m2sm_map{ SubMeshUtils::BuildFaceMap(m, tf_sm, tf_sm.GetParentElementIDMap()) };
	auto sf_m2sm_map{ SubMeshUtils::BuildFaceMap(m, sf_sm, sf_sm.GetParentElementIDMap()) };
	auto f2bdr_map{ m.GetFaceToBdrElMap() };
	for (int i = 0; i < m.GetNBE(); i++) {
		if (m.GetBdrAttribute(i) == 301) {
			tf_sm.SetBdrAttribute(tf_m2sm_map.Find(f2bdr_map.Find(i)), 301);
			sf_sm.SetBdrAttribute(sf_m2sm_map.Find(f2bdr_map.Find(i)), 301);
		}
	}

	std::vector<int> tf_bdr_ids{ {0,1,2,3} };
	std::vector<int> sf_bdr_ids{ {1,6,8,15} };

	for (int i = 0; i < 4; i++) {
		EXPECT_TRUE(tf_sm.GetBdrAttribute(tf_bdr_ids[i]) == 301);
		EXPECT_TRUE(sf_sm.GetBdrAttribute(sf_bdr_ids[i]) == 301);
	}
}

TEST_F(GeometryTest, original_attributes_are_restored_2D)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square3x3marked.mesh").c_str(), 1, 0, true) };

	setTFSFAttributesForSubMeshing(m);

	Array<int> tf_att(1); tf_att[0] = 1000;
	Array<int> sf_att(1); sf_att[0] = 2000;
	auto tf_sm{ SubMesh::CreateFromDomain(m, tf_att) };
	auto sf_sm{ SubMesh::CreateFromDomain(m, sf_att) };

	restoreOriginalAttributes(m);

	for (int e = 0; e < m.GetNE(); e++) {
		EXPECT_TRUE(m.GetElement(e)->GetAttribute() == 1);
	}
}



