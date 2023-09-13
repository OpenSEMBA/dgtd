#include <gtest/gtest.h>

#include <mfem.hpp>

#include "TestUtils.h"  

using namespace mfem;

using NodeId = int;
using FaceId = int;
using ElementId = int;
using Orientation = int;
using Attribute = int;
using BdrId = int;
using IsInterior = bool;
using ElNo2Att = std::pair<ElementId, Attribute>;
using TwoElems = std::pair<ElementId, ElementId>;
using FaceToAtt = std::map<FaceId, Attribute>;

double linearFunction(const Vector& pos)
{
	double normalizedPos;
	double leftBoundary = 0.0, rightBoundary = 1.0;
	double length = rightBoundary - leftBoundary;
	normalizedPos = (pos[0] - leftBoundary) / length;

	return 2 * normalizedPos;
}

std::map<FaceId, Orientation> mapFaceToOrientationOuterBoundary(const Mesh& mesh, const size_t numberOfElements) {
	
	/*This method would require from another method that would identify elements with a specific tag,
	then the numberOfElements argument would be removed and the variable would be substituted by
	the # of tagged elements.*/
	
	std::map<FaceId, Orientation> res;
	for (auto i{ 0 }; i < numberOfElements; ++i) {
		Array<FaceId> faces;
		Array<Orientation> orientations;
		mesh.GetElementEdges(i, faces, orientations);

		assert(faces.Size() == orientations.Size());
		for (auto f{ 0 }; f < faces.Size(); ++f) {
			auto it{ res.find(faces[f]) };
			if (it == res.end()) {
				res[faces[f]] = orientations[f];
			}
			else {
				res.erase(it);
			}
		}
	}
	return res;
}

class MeshTest : public ::testing::Test {
protected:
	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian1D(1);
	}

	Mesh mesh_;

	static std::string getFilename(const std::string fn)
	{
		return "./testData/" + fn;
	}

	std::vector<int> mapQuadElementTopLeftVertex(const Mesh& mesh)
	{
		std::vector<int> res;
		for (int i = 0; i < mesh.GetNE(); i++) {
			Array<int> meshArrayElement;
			mesh.GetElementVertices(i, meshArrayElement);
			res.push_back(meshArrayElement[0]);
		}
		return res;
	}

	Mesh makeTwoAttributeCartesianMesh1D(const int& refTimes = 0)
	{
		Mesh res = Mesh::MakeCartesian1D(2);
		res.SetAttribute(0, 1);
		res.SetAttribute(1, 2);

		for (int i = 0; i < refTimes; i++) {
			res.UniformRefinement();
		}

		return res;
	}


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

TEST_F(MeshTest, TwoAttributeMesh)
{
	/*The purpose of this test is to check the makeTwoAttributeCartesianMesh1D(const int& refTimes)
	function.

	First, an integer is declared for the number of times we wish to refine the mesh, then a mesh is
	constructed with two elements, left and right hand sides, setting the following attributes.

	|------LHS------|------RHS------|

	|##ATTRIBUTE 1##|##ATTRIBUTE 2##|

	Once the mesh is refined, it is returned, then we compare if the expected number of elements is
	true for the actual elements in the mesh.

	Then, we consider how the mesh will perform its uniform refinement, and we declare that the
	LHS elements with Attribute one will be Even index elements (starting at 0), and the RHS
	elements with Attribute 2 will be Uneven index elements (starting at 1).*/

	const int refTimes = 3;
	Mesh mesh = makeTwoAttributeCartesianMesh1D(refTimes);

	EXPECT_EQ(pow(2, refTimes + 1), mesh.GetNE());
	for (int i = 0; i < mesh.GetNE(); i++) {
		if (i % 2 == 0) {
			EXPECT_EQ(1, mesh.GetAttribute(i));
		}
		else {
			EXPECT_EQ(2, mesh.GetAttribute(i));
		}
	}
}

TEST_F(MeshTest, MeshDimensions)
{
	/*This test ensures that the number of elements of any 2D Cartesian
	mesh is equal to the product of the horizontal and vertical segments

	Dimensional parameters are declared which are then used to create
	a mesh object, then the test comparison is made.*/


	int nx = 8; int ny = 8; bool generateEdges = true;
	Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, generateEdges);

	EXPECT_EQ(nx * ny, mesh.GetNE());

}

TEST_F(MeshTest, ElementVolume_2D_Triangle)
{

	Mesh mesh{ Mesh::MakeCartesian2D(1,1,Element::TRIANGLE)};
	for (int i = 0; i < mesh.GetNE(); ++i) {
		ASSERT_EQ(0.5, mesh.GetElementVolume(i));
	}

}

TEST_F(MeshTest, ElementVolume_3D_Tetra)
{
	Mesh mesh{ Mesh::MakeCartesian3D(1,1,1,Element::TETRAHEDRON) };
	for (int i = 0; i < mesh.GetNE(); ++i) {
		ASSERT_EQ(1.0 / 6.0, mesh.GetElementVolume(i));
	}
}

TEST_F(MeshTest, ElementFaceSurface_3D_Tetra)
{
	Mesh mesh{ Mesh::MakeCartesian3D(1,1,1,Element::TETRAHEDRON) };
	Vector surfaces(mesh.GetNFaces());
	for (int i = 0; i < mesh.GetNE(); ++i) {
		Array<int> faces, ori;
		mesh.GetElementFaces(i, faces, ori);
		for (int j = 0; j < faces.Size(); ++j) {
			auto face{ mesh.GetFace(faces[j]) };
			auto T{ mesh.GetFaceElementTransformations(faces[j]) };
			IntegrationRule ir;
			double res = 0.0;
			for (int p = 0; p < ir.GetNPoints(); p++)
			{
				const IntegrationPoint& ip = ir.IntPoint(p);
				T->SetAllIntPoints(&ip);
				res += ip.weight * T->Weight();
			}
			surfaces(faces[j]) = res;
		}
	}
	//WIP
}

TEST_F(MeshTest, ElementPerimeter_2D_Triangle)
{
	Mesh mesh{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };
	Vector per(mesh.GetNE());
	per = 0.0;

	for (int i = 0; i < mesh.GetNE(); ++i) {
		auto el{ mesh.GetElement(i) };
		for (int f = 0; f < mesh.GetElement(i)->GetNEdges(); ++f) {
			ElementTransformation* T{ mesh.GetFaceTransformation(f) };
			const IntegrationRule& ir = IntRules.Get(T->GetGeometryType(), T->OrderJ());
			for (int p = 0; p < ir.GetNPoints(); p++)
			{
				const IntegrationPoint& ip = ir.IntPoint(p);
				per(i) += ip.weight * T->Weight();
			}
		}
	}

	double tol{ 1e-8 };
	for (int i = 0; i < mesh.GetNE(); ++i) {
		ASSERT_NEAR(1.0 + 1.0 + sqrt(pow(1.0,2.0) + pow(1.0,2.0)), per(i), tol);
	}

}

TEST_F(MeshTest, ElementVolumeThroughPerimeter_2D_Triangle)
{
	Mesh mesh{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };

	auto NV{ mesh.GetNV() };
	Vector vertCoord(NV);
	mesh.GetVertices(vertCoord);
	Vector vx(NV), vy(NV);
	for (int i = 0; i < NV; ++i) {
		vx(i) = vertCoord(i);
		vy(i) = vertCoord(i + NV);
	}

	Vector area(mesh.GetNE());
	for (int it = 0; it < mesh.GetNE(); ++it) {
		auto el{ mesh.GetElement(it) };
		Array<int> EV(el->GetNVertices());
		el->GetVertices(EV);
		Vector len(EV.Size());
		len = 0.0;
		double sper{ 0.0 };
		for (int i = 0; i < EV.Size(); ++i) {
			int j = (i + 1) % 3;
			len(i) += sqrt(pow(vx(EV[i]) - vx(EV[j]), 2.0) + pow(vy(EV[i]) - vy(EV[j]), 2.0));
			sper += len(i);
		}
		sper /= 2.0;
		area(it) = sqrt(sper);
		for (int i = 0; i < EV.Size(); ++i) {
			area(it) *= sqrt(sper - len(i));
		}
	}

	double tol = 1e-8;
	for (int i = 0; i < mesh.GetNE(); ++i) {
		ASSERT_NEAR(mesh.GetElementVolume(i), area(i),tol);
	}

}

TEST_F(MeshTest, EstimatedDTScale_2D_Triangle) 
{
	Mesh mesh{ Mesh::MakeCartesian2D(1,1,Element::Type::TRIANGLE) };

	Vector area(mesh.GetNE()), dtscale(mesh.GetNE());
	for (int it = 0; it < mesh.GetNE(); ++it) {
		auto el{ mesh.GetElement(it) };
		Vector sper(mesh.GetNumFaces());
		sper = 0.0;
		for (int f = 0; f < mesh.GetElement(it)->GetNEdges(); ++f) {
			ElementTransformation* T{ mesh.GetFaceTransformation(f) };
			const IntegrationRule& ir = IntRules.Get(T->GetGeometryType(), T->OrderJ());
			for (int p = 0; p < ir.GetNPoints(); p++)
			{
				const IntegrationPoint& ip = ir.IntPoint(p);
				sper(it) += ip.weight * T->Weight();
			}
		}
		sper /= 2.0;
		area(it) = mesh.GetElementVolume(it);
		dtscale(it) = area(it) / sper(it);
	}

	double tol = 1e-8;
	for (int i = 0; i < mesh.GetNE(); ++i) {
		ASSERT_NEAR(0.292893218813452, dtscale(i), tol);
	}

}

TEST_F(MeshTest, DataValueOutsideNodesForOneElementMeshes)
{
	/* The purpose of this test is to ensure we can extract data from a GridFunction,
	even if the point we're trying to obtain it at is not necessarily a DoF or node.
	
	First, the basic process to declare and initialise a FiniteElementSpace is done,
	this means variables such as dimension, order, creating a mesh, a FEC and a finally,
	the FES.
	
	A GridFunction is then created and assigned the FES. A function is projected in the
	GridFunction, which is a linear function with a slope of 2.
	
	Lastly, an IntegrationPoint is constructed, which we will use to obtain the values
	from the GridFunction at any point we want. As the slope of the line is 2, we expect
	the values to be 2 times the xVal.*/
	
	Mesh mesh = Mesh::MakeCartesian1D(1);
	auto fecDG = new DG_FECollection(1, 1, BasisType::GaussLobatto);
	auto* fesDG = new FiniteElementSpace(&mesh, fecDG);

	GridFunction solution{ fesDG };
	FunctionCoefficient fc{ linearFunction };
	solution.ProjectCoefficient(fc);
	IntegrationPoint integPoint;
	for (double xVal = 0.0; xVal <= 1; xVal = xVal + 0.1) {
		integPoint.Set(xVal, 0.0, 0.0, 0.0);
		double interpolatedPoint = solution.GetValue(0, integPoint);
		EXPECT_NEAR(xVal * 2, interpolatedPoint,1e-10);
	}
}

TEST_F(MeshTest, MeshElementVertices)
{
	/*This test was created to understand the process of mesh creation
	and assignation of vertex index to elements.

	First, dimensional variables are declared, which are then used
	to create a new mesh object.

	Then firstElementVerticesVector and lastElementVerticesVector are
	initialized and assigned values manually, with the values we expect
	these elements will have. We also retrieve the vertices for the first
	and last element of the mesh, and then store them inside vectors.

	Lastly, we compare that the vertices we retrieved from the mesh are
	equal to those we presumed at the start.*/

	int nx = 8; int ny = 8; bool generateEdges = true;
	Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, generateEdges);

	std::vector<int> firstElementVerticesVector = { 0, 1, nx + 2, nx + 1 };
	std::vector<int> lastElementVerticesVector = { nx - 1, nx, nx * 2 + 1, nx * 2 };
	Array<int> meshArrayFirstElement;
	Array<int> meshArrayLastElement;

	mesh.GetElementVertices(0, meshArrayFirstElement);
	mesh.GetElementVertices(nx * ny - 1, meshArrayLastElement);

	std::vector<int> vectorFirstElement(meshArrayFirstElement.begin(), meshArrayFirstElement.end());
	std::vector<int> vectorLastElement(meshArrayLastElement.begin(), meshArrayLastElement.end());

	EXPECT_EQ(firstElementVerticesVector, vectorFirstElement);
	EXPECT_EQ(lastElementVerticesVector, vectorLastElement);

}

TEST_F(MeshTest, MapMeshElementAndVertex)
{

	/* This test was created with the aim to understand the mapping and ordering process
	of a mesh in a more visual way. It uses the mapQuadElementTopLeftVertex() function
	which, for a Quadrilateral Element, it extracts its top left vertex, which allows for a nigh
	full mapping of the mesh.

	First, dimensional variables are declared and a mesh is constructed.

	Then, the mapQuadElementTopLeftVertex extracts the top left vertex of each element and stores them
	in an integer vector.

	Lastly, we compare that the first mapped vertex is the first created vertex in the mesh 0,
	the top left vertex for the uppermost, rightmost element is equal to the last element's index - 1
	(due to how mesh mapping works), and the size of the mapped vertices vector is equal to the number
	of elements in the mesh - 1 (as it starts with index 0).*/

	int nx = 5; int ny = 5; bool generateEdges = true;
	Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, generateEdges);

	std::vector<int> mapped = mapQuadElementTopLeftVertex(mesh);

	EXPECT_EQ(0, mapped[0]);
	EXPECT_EQ(nx - 1, mapped[mapped.size() - 1]);
	EXPECT_EQ(nx * ny - 1, mapped.size() - 1);

}

TEST_F(MeshTest, MeshDataFileRead)
{
	ASSERT_NO_THROW(Mesh::LoadFromFile((mfemMeshesFolder() + "twotriang.mesh").c_str(), 1, 0));
}

TEST_F(MeshTest, BoundaryWithoutInteriorFace)
{
	const auto numberOfElements{ 2 };
	auto mesh{ Mesh::MakeCartesian2D(numberOfElements, 1, Element::Type::QUADRILATERAL, false, 2.0) };

	auto facesToOrient = mapFaceToOrientationOuterBoundary(mesh, numberOfElements);

	std::map<int, int> expected{{0,1},{2,-1},{3,-1},{4,1},{5,1},{6,-1}};

	EXPECT_EQ(expected, facesToOrient);
}

TEST_F(MeshTest, SubMeshingAttributes_1D)
{
	auto mesh{ Mesh::MakeCartesian1D(2) };

	const auto att{ 3 };
	Array<int> subdomain_attributes(1);
	subdomain_attributes[0] = att;

	for (int i = 0; i < mesh.GetNE(); ++i) {
		const auto preAtt{ mesh.GetAttribute(i) };
		mesh.SetAttribute(i, att);
		auto submesh{ SubMesh::CreateFromDomain(mesh, subdomain_attributes) };
		mesh.SetAttribute(i, preAtt);
		EXPECT_EQ(1,mesh.GetAttribute(i));
		EXPECT_EQ(3,submesh.GetAttribute(0));
	}

}

TEST_F(MeshTest, SubMeshingBdrAttributes_1D)
{
	auto mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "line_for_submesh.mesh").c_str(),1, 0) };

	Array<int> subdomain_attribute_1(1); subdomain_attribute_1[0] = 1;
	Array<int> subdomain_attribute_2(1); subdomain_attribute_2[0] = 2;
	Array<int> subdomain_attribute_3(1); subdomain_attribute_3[0] = 3;

	auto submesh_att_1{ SubMesh::CreateFromDomain(mesh, subdomain_attribute_1) };
	auto submesh_att_2{ SubMesh::CreateFromDomain(mesh, subdomain_attribute_2) };
	auto submesh_att_3{ SubMesh::CreateFromDomain(mesh, subdomain_attribute_3) };

	submesh_att_1.AddBdrPoint(0, 4);
	submesh_att_3.AddBdrPoint(1, 5);

	EXPECT_EQ(mesh.GetAttribute(0), submesh_att_1.GetAttribute(0));
	EXPECT_EQ(mesh.GetAttribute(1), submesh_att_2.GetAttribute(0));
	EXPECT_EQ(mesh.GetAttribute(2), submesh_att_3.GetAttribute(0));
	EXPECT_EQ(mesh.GetBdrAttribute(0), submesh_att_1.GetBdrAttribute(0));
	EXPECT_EQ(mesh.GetBdrAttribute(1), submesh_att_3.GetBdrAttribute(0));

}

TEST_F(MeshTest, MeshIdentifyBoundaryVertex)
{
	auto mesh{ Mesh::MakeCartesian1D(2) };

	Vector bdrIndex(mesh.GetNBE());
	for (int i = 0; i < bdrIndex.Size(); ++i) {
		bdrIndex(i) = int(mesh.GetBdrElementEdgeIndex(i));
	}
	
	Vector exp({ 0,2 });
	for (int i = 0; i < bdrIndex.Size(); ++i) {
		EXPECT_EQ(exp[0], bdrIndex[0]);
	}

}

TEST_F(MeshTest, SubMeshAssignBdrFromParentMesh_1D)
{
	auto mesh{ Mesh::MakeCartesian1D(2) };
	mesh.SetBdrAttribute(0, 2);	mesh.SetBdrAttribute(1, 3);

	using ParentId = int;
	using BdrAtt = int;

	Vector bdrIndex(mesh.GetNBE());
	for (int i = 0; i < bdrIndex.Size(); ++i) {
		bdrIndex(i) = int(mesh.GetBdrElementEdgeIndex(i));
	}

	for (int e = 0; e < mesh.GetNE(); ++e) {
		std::map<ParentId, BdrAtt> PIdToBdrAtt;
		auto el{ mesh.GetElement(e) };
		Array<int> ver(el->GetNVertices());
		Array<int> bdrAtt(el->GetNVertices());
		Array<bool> isBdr(el->GetNVertices());
		el->GetVertices(ver);
		for (int v = 0; v < ver.Size(); ++v) {
			mesh.FaceIsInterior(ver[v]) == false ? isBdr[v] = true : isBdr[v] = false;
			mesh.FaceIsInterior(ver[v]) == false ? bdrAtt[v] = 0   : bdrAtt[v] = 3;
		}
	}

}

TEST_F(MeshTest, SubMeshingOnX)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "FiveQuadsOnX_TFSF.mesh").c_str(), 1, 0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLobatto) };
	auto fes_m{ FiniteElementSpace(&m,&fec) };

	Array<int> domain_atts({ 201, 202 });
	auto sm{ SubMesh::CreateFromDomain(m,domain_atts) };
	Array<int> att201, att202;
	for (int e = 0; e < sm.GetNE(); ++e) {
		if (sm.GetAttribute(e) == 201) {
			att201.Append(e);
		}
		if (sm.GetAttribute(e) == 202) {
			att202.Append(e);
		}
	}
	MFEM_ASSERT(att201 == Array<int>({0, 3}), "att201 failure.");
	MFEM_ASSERT(att202 == Array<int>({1, 2}), "att202 failure.");
	
	auto fes_sm{ FiniteElementSpace(&sm,&fec) };
	
	GridFunction gf_m (&fes_m);
	GridFunction gf_sm(&fes_sm);

	gf_sm[1]  =  5.0;
	gf_sm[13] = -5.0;

	sm.Transfer(gf_sm, gf_m);

	auto idmap{ sm.GetParentElementIDMap() };
	FaceToAtt parentbdr2a;
	for (int bdr = 0; bdr < m.GetNBE(); ++bdr) {
		parentbdr2a.emplace(bdr, m.GetBdrAttribute(bdr));
	}
	FaceToAtt subbdr2a;
	for (int bdr = 0; bdr < sm.GetNBE(); ++bdr) {
		subbdr2a.emplace(bdr, sm.GetBdrAttribute(bdr));
	}

	auto parentf2bdrEl{ m.GetFaceToBdrElMap() };
	auto subf2bdrEl{ sm.GetFaceToBdrElMap() };

	//Overrider of Boundaries, Squasher of Ones, Reclassifier of Sources

	for (int e = 0; e < subf2bdrEl.Size(); ++e) { //Not all non-bdrs have to be internal or external bdrs, needs rework.
		if (subf2bdrEl[e] != -1) { continue; }
		Array<int> v;
		sm.GetEdgeVertices(e, v);
		sm.AddBdrSegment(v[0], v[1]);
	}

	auto m2smFaceMap{ mfem::SubMeshUtils::BuildFaceMap(m, sm, idmap) };
	for (int f = 0; f < parentf2bdrEl.Size(); ++f) {
		if (parentf2bdrEl[f] == -1) { continue; }
		auto attToBe{ parentbdr2a[parentf2bdrEl[f]] };

	}


}

TEST_F(MeshTest, SubMeshingFromDomain)
{
	auto m{ Mesh::LoadFromFile((gmshMeshesFolder() + "TotalFieldScatteredFieldQuads.msh").c_str(), 1, 0)};
	
	Array<int> domain_atts(2); domain_atts[0] = 201; domain_atts[1] = 202;


	std::map<ElementId, FaceToAtt> map_pair;
	std::map<BdrId, Attribute> bdr2att_map;
	std::map<BdrId, TwoElems> bdr2el_map;
	std::map<BdrId, IsInterior> bdr2int_map;
	bdr2int_map.emplace(2, false); bdr2int_map.emplace(301, true);

	for (int b = 0; b < m.GetNBE(); ++b) {
		bdr2att_map.emplace(b,m.GetBdrAttribute(b));
		int el, el2, info, info2;
		if (bdr2int_map[m.GetBdrAttribute(b)] == true) {
			auto facetrans = m.GetFaceElementTransformations(b, 31);
			el = facetrans->Elem1No;
			el2 = facetrans->Elem2No;
		}
		else {
			auto facetrans = m.GetBdrFaceTransformations(b);
			el = facetrans->Elem1No;
			el2 = -1;
		}
		TwoElems twoelems(el,el2);
		bdr2el_map.emplace(b,twoelems);
	}

	auto faceneighs{ m.FindFaceNeighbors(2) };	

	for (int e = 0; e < m.GetNE(); ++e) {
		for (int i = 0; i < domain_atts.Size(); ++i) {
			if (m.GetElement(e)->GetAttribute() == domain_atts[i]) {
				Array<int> faces, ori;
				for (int f = 0; f < m.GetElement(e)->GetNEdges(); ++f) {
					m.GetElementFaces(e, faces, ori);
				}

			}
		}
	}

	auto sm_dom{ SubMesh::CreateFromDomain(m,domain_atts) };
	Array<int> bdr_marker(2); bdr_marker[0] = 2; bdr_marker[1] = 301;
	auto sm_bdr{ SubMesh::CreateFromBoundary(m,bdr_marker) };


}

TEST_F(MeshTest, ExtendedFindNeighboursMethod2D)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "square3x3.mesh").c_str(), 1, 0) };
	
	MFEM_ASSERT(m.FindFaceNeighbors(0) == Array<int>({0, 1, 3})      ,"Elem 0 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(1) == Array<int>({0, 1, 2, 4})   ,"Elem 1 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(2) == Array<int>({1, 2, 5})      ,"Elem 2 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(3) == Array<int>({0, 3, 4, 6})   ,"Elem 3 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(4) == Array<int>({1, 3, 4, 5, 7}),"Elem 4 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(5) == Array<int>({2, 4, 5, 8})   ,"Elem 5 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(6) == Array<int>({3, 6, 7})      ,"Elem 6 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(7) == Array<int>({4, 6, 7, 8})   ,"Elem 7 fails.");
	MFEM_ASSERT(m.FindFaceNeighbors(8) == Array<int>({5, 7, 8})      ,"Elem 8 fails.");
}

TEST_F(MeshTest, InteriorBoundary)
{
	auto mesh{ Mesh::LoadFromFile((mfemMeshesFolder() + "line.mesh").c_str(), 1, 0) };

	EXPECT_EQ(2, mesh.GetBdrAttribute(0));
	EXPECT_EQ(2, mesh.GetBdrAttribute(1));
	EXPECT_EQ(3, mesh.GetBdrAttribute(2));

	mesh.AddBdrPoint(2, 4);
	EXPECT_EQ(4, mesh.GetBdrAttribute(3));
}

TEST_F(MeshTest, GetElementSize_1D)
{
	auto m{ Mesh::MakeCartesian1D(10, 1.0) };

	for (auto e{ 0 }; e < m.GetNE(); ++e) {
		EXPECT_NEAR(0.1, m.GetElementSize(e), 1e-8);
	}

}

TEST_F(MeshTest, GetElementFaceOrientation)
{
	/*    2___-1__3
		  | 	  |
		  -1	  1
		  |		  |
	      0___1___1
	*/


	auto m{ Mesh::MakeCartesian2D(1,1,Element::QUADRILATERAL,true) };

	Array<int> faces(m.GetElement(0)->GetNEdges()), ori(m.GetElement(0)->GetNEdges());
	m.GetElementEdges(0, faces, ori);
	Array<int> vert0(2), vert1(2), vert2(2), vert3(2);
	m.GetEdgeVertices(0,vert0);
	m.GetEdgeVertices(1,vert1);
	m.GetEdgeVertices(2,vert2);
	m.GetEdgeVertices(3,vert3);
}

TEST_F(MeshTest, FaceElementSurface_2D_Tri)
{
	auto m{ Mesh::MakeCartesian2D(1,1,Element::TRIANGLE) };
	Vector surface(m.GetNumFaces());
	surface = 0.0;
	for (int f = 0; f < m.GetNumFaces(); ++f) {
		ElementTransformation* T{ m.GetFaceTransformation(f) };
		const IntegrationRule& ir = IntRules.Get(T->GetGeometryType(), T->OrderJ());
		for (int p = 0; p < ir.GetNPoints(); p++)
		{
			const IntegrationPoint& ip = ir.IntPoint(p);
			surface(f) += ip.weight * T->Weight();
		}
	}
	for (int e = 0; e < m.GetNE(); ++e) {
		Array<int> edges(m.GetElement(e)->GetNEdges()), ori(m.GetElement(e)->GetNEdges());
		m.GetElementEdges(e, edges, ori);
		double tol = 1e-8;
		ASSERT_NEAR(1.0 + 1.0 + sqrt(2), surface(edges[0]) + surface(edges[1]) + surface(edges[2]), tol);
	}
}

TEST_F(MeshTest, ElementEdges_2D_Tri)
{

/* For a triangular mesh with elements and edges such as
* 
*                  ____1_____
*                 |        / |
*                 |	 0   /   |
*                 2    0     4
*                 |	 /    1  |
*                 |/___3_____|
*
* we make a test where we obtain the edges of each element
* to understand expected edge ordering and surface value.*/

	auto m{ Mesh::MakeCartesian2D(1,1,Element::TRIANGLE) };
	for (int e = 0; e < m.GetNE(); ++e) {
		Array<int> edges(m.GetElement(e)->GetNEdges()), ori(m.GetElement(e)->GetNEdges());
		m.GetElementEdges(e, edges, ori);
		if (e == 0) {

			ASSERT_EQ(Array<int>({0, 1, 2 }), edges);
		}
		else {
			ASSERT_EQ(Array<int>({0, 3, 4 }), edges);
		}
	}
}

TEST_F(MeshTest, SubMeshingAttributes_2D)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "four_quads_for_submeshing.mesh").c_str(),1,0) };

	Array<int> att_1(1); att_1[0] = 1;
	Array<int> att_2(1); att_2[0] = 2;
	Array<int> att_3(1); att_3[0] = 3;
	Array<int> att_4(1); att_4[0] = 4;

	auto submesh_att_1{ SubMesh::CreateFromDomain(m,att_1)};
	auto submesh_att_2{ SubMesh::CreateFromDomain(m,att_2)};
	auto submesh_att_3{ SubMesh::CreateFromDomain(m,att_3)};
	auto submesh_att_4{ SubMesh::CreateFromDomain(m,att_4)};

	EXPECT_EQ(m.GetAttribute(0), submesh_att_1.GetAttribute(0));
	EXPECT_EQ(m.GetAttribute(1), submesh_att_2.GetAttribute(0));
	EXPECT_EQ(m.GetAttribute(2), submesh_att_3.GetAttribute(0));
	EXPECT_EQ(m.GetAttribute(3), submesh_att_4.GetAttribute(0));
	EXPECT_EQ(1, submesh_att_1.GetAttribute(0));
	EXPECT_EQ(2, submesh_att_2.GetAttribute(0));
	EXPECT_EQ(3, submesh_att_3.GetAttribute(0));
	EXPECT_EQ(4, submesh_att_4.GetAttribute(0));

}

TEST_F(MeshTest, SubMeshingBdrAttributes_2D)
{
	auto m{ Mesh::LoadFromFile((mfemMeshesFolder() + "four_quads_for_submeshing.mesh").c_str(),1,0) };

	Array<int> att_1(1); att_1[0] = 1;
	Array<int> att_2(1); att_2[0] = 2;
	Array<int> att_3(1); att_3[0] = 3;
	Array<int> att_4(1); att_4[0] = 4;

	auto submesh_bdratt_0_1{ SubMesh::CreateFromDomain(m,att_1) };
	auto submesh_bdratt_2_3{ SubMesh::CreateFromDomain(m,att_2) };
	auto submesh_bdratt_4_5{ SubMesh::CreateFromDomain(m,att_3) };
	auto submesh_bdratt_6_7{ SubMesh::CreateFromDomain(m,att_4) };

	EXPECT_EQ(m.GetBdrAttribute(0), submesh_bdratt_0_1.GetBdrAttribute(0));
	EXPECT_EQ(m.GetBdrAttribute(1), submesh_bdratt_0_1.GetBdrAttribute(1));
	EXPECT_EQ(m.GetBdrAttribute(2), submesh_bdratt_2_3.GetBdrAttribute(0));
	EXPECT_EQ(m.GetBdrAttribute(3), submesh_bdratt_2_3.GetBdrAttribute(1));
	EXPECT_EQ(m.GetBdrAttribute(4), submesh_bdratt_4_5.GetBdrAttribute(0));
	EXPECT_EQ(m.GetBdrAttribute(5), submesh_bdratt_4_5.GetBdrAttribute(1));
	EXPECT_EQ(m.GetBdrAttribute(6), submesh_bdratt_6_7.GetBdrAttribute(0));
	EXPECT_EQ(m.GetBdrAttribute(7), submesh_bdratt_6_7.GetBdrAttribute(1));
	EXPECT_EQ(5, submesh_bdratt_0_1.GetBdrAttribute(0));
	EXPECT_EQ(5, submesh_bdratt_0_1.GetBdrAttribute(1));
	EXPECT_EQ(6, submesh_bdratt_2_3.GetBdrAttribute(0));
	EXPECT_EQ(6, submesh_bdratt_2_3.GetBdrAttribute(1));
	EXPECT_EQ(7, submesh_bdratt_4_5.GetBdrAttribute(0));
	EXPECT_EQ(7, submesh_bdratt_4_5.GetBdrAttribute(1));
	EXPECT_EQ(8, submesh_bdratt_6_7.GetBdrAttribute(0));
	EXPECT_EQ(8, submesh_bdratt_6_7.GetBdrAttribute(1));
	
}


TEST_F(MeshTest, marking_element_att_through_boundary_2D)
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

TEST_F(MeshTest, marking_element_to_face_pairs_for_submeshing_2D)
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

TEST_F(MeshTest, transfer_parent_bdr_att_to_child)
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

TEST_F(MeshTest, original_attributes_are_restored_2D)
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
