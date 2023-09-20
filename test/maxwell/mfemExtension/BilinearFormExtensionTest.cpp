#include <gtest/gtest.h>

#include "mfemExtension/BilinearIntegrators.h"
#include "components/Types.h"
#include <TestUtils.h>

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

using FaceId = int;
using ElementId = int; 
using Attribute = int; 
using ElNo2Att = std::pair<ElementId, Attribute>;


class BilinearFormExtensionTest : public ::testing::Test 
{
protected:

	void SetUp() override
	{
		mesh_ = Mesh::MakeCartesian1D(1);
		fec_ = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES1D(
		const int order,
		const int elements = 1,
		const double length = 1.0,
		const decltype(BasisType::GaussLobatto) basis = BasisType::GaussLobatto
	)
	{
		mesh_ = Mesh::MakeCartesian1D(elements, length);
		fec_ = std::make_unique<DG_FECollection>(order, 1, basis);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

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

TEST_F(BilinearFormExtensionTest, buildBilinearFormFromSubMeshes)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square3x3marked.mesh").c_str(), 1, 0, true) };

	setTFSFAttributesForSubMeshing(m);

	Array<int> tf_att(1); tf_att[0] = 1000;
	Array<int> sf_att(1); sf_att[0] = 2000;
	auto tf_sm{ SubMesh::CreateFromDomain(m, tf_att) };
	auto sf_sm{ SubMesh::CreateFromDomain(m, sf_att) };
	auto tf_m2sm_map{ SubMeshUtils::BuildFaceMap(m, tf_sm, tf_sm.GetParentElementIDMap()) };
	auto sf_m2sm_map{ SubMeshUtils::BuildFaceMap(m, sf_sm, sf_sm.GetParentElementIDMap()) };
	auto f2bdr_map{ m.GetFaceToBdrElMap() };
	for (int i = 0; i < m.GetNBE(); i++) {
		if (m.GetBdrAttribute(i) == 301) {
			tf_sm.SetBdrAttribute(tf_m2sm_map.Find(f2bdr_map.Find(i)), 301);
			sf_sm.SetBdrAttribute(sf_m2sm_map.Find(f2bdr_map.Find(i)), 301);
		}
	}
	tf_sm.FinalizeMesh();
	sf_sm.FinalizeMesh();

	auto fec{ L2_FECollection{1, 2, BasisType::GaussLobatto }};
	auto tf_fes{ FiniteElementSpace{&tf_sm, &fec} };

	BilinearForm tf_bf(&tf_fes);
	Array<int> bdr_marker(301);
	bdr_marker = 0; bdr_marker[300] = 1;
	tf_bf.AddBdrFaceIntegrator(new TotalFieldScatteredFieldIntegrator(1.0), bdr_marker);
	tf_bf.Assemble();
	tf_bf.Finalize();
	tf_bf.SpMat().ToDenseMatrix()->Print(std::cout);
	std::cout << std::flush;

}

//TEST_F(BilinearFormExtensionTest, checkInteriorBoundaryFaceIntegrator)
//{
//	setFES1D(1,4,4.0);
//
//	auto intBdrAttr{ 5 };
//	mesh_.AddBdrPoint(2, intBdrAttr);
//	mesh_.FinalizeMesh();
//
//	BilinearFormIBFI totalFieldFlux{ fes_.get() };
//	Array<int> intBdrMarker{ mesh_.bdr_attributes.Max() };
//	intBdrMarker = 0;
//	intBdrMarker[intBdrAttr - 1] = 1;
//	std::vector<VectorConstantCoefficient> n = {VectorConstantCoefficient(Vector({1.0}))};
//	totalFieldFlux.AddInteriorBoundaryFaceIntegrator(
//		new DGTraceIntegrator{n[0],0.0, 1.0},
//		intBdrMarker
//	);
//	totalFieldFlux.Assemble();
//	totalFieldFlux.Finalize();
//
//	GridFunction f{ fes_.get() }, exc{ fes_.get() };
//	f = 0.0; 
//	exc[3] = 6.0;
//	exc[4] = 5.0;
//	
//	totalFieldFlux.Mult(exc, f);
//
//	EXPECT_EQ( 0.0, f[0]);
//	EXPECT_EQ( 0.0, f[7]);
//	EXPECT_EQ( 1.0, f[3]);
//	EXPECT_EQ(-1.0, f[4]);
//}

//TEST_F(BilinearFormExtensionTest, compareBaseAndDerivedBilinearForms)
//{
//	setFES1D(1, 4, 4.0);
//
//	auto intBdrAttr{ 5 };
//	mesh_.AddBdrPoint(2, intBdrAttr);
//	mesh_.FinalizeMesh();
//
//	Array<int> intBdrMarker{ mesh_.bdr_attributes.Max() };
//	intBdrMarker = 0;
//	intBdrMarker[intBdrAttr - 1] = 1;
//	std::vector<VectorConstantCoefficient> n = { VectorConstantCoefficient(Vector({1.0})) };
//
//	BilinearForm baseBilinearForm(fes_.get());
//	BilinearFormIBFI derivedBilinearForm(fes_.get());
//	baseBilinearForm.AddBdrFaceIntegrator(
//		new DGTraceIntegrator{ n[0],0.0, 1.0 },
//		intBdrMarker
//	);
//	derivedBilinearForm.AddInteriorBoundaryFaceIntegrator(
//		new DGTraceIntegrator{ n[0],0.0, 1.0 },
//		intBdrMarker
//	);
//	baseBilinearForm.Assemble();
//	baseBilinearForm.Finalize();
//	derivedBilinearForm.Assemble();
//	derivedBilinearForm.Finalize();	
//
//	GridFunction fbase{ fes_.get() }, fderived{ fes_.get() },exc{ fes_.get() };
//	exc[3] = 6.0;
//	exc[4] = 5.0;
//	baseBilinearForm.Mult(exc, fbase);
//	derivedBilinearForm.Mult(exc, fderived);
//
//	EXPECT_EQ( 0.0, fbase[3]);
//	EXPECT_EQ( 0.0, fbase[4]);
//	EXPECT_EQ( 1.0, fderived[3]);
//	EXPECT_EQ(-1.0, fderived[4]);
//	
//}


