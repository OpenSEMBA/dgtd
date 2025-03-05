#include <gtest/gtest.h>

#include <mfem.hpp>

#include "TestUtils.h"  
#include "components/SubMesher.h"
#include "components/Types.h"

#include "evolution/HesthavenEvolutionMethods.h"
#include <evolution/HesthavenEvolution.h>

using namespace mfem;
using namespace maxwell;

using NodeId = int;
using FaceId = int;
using ElementId = int;
using Orientation = int;
using Attribute = int;
using BdrId = int;
using IsInterior = bool;
using IsTF = bool;
static int NotFound{ -1 };
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
	std::unique_ptr<DG_FECollection> fecDG = std::make_unique<DG_FECollection>(1, 1, BasisType::GaussLobatto);
	auto fesDG = FiniteElementSpace(&mesh, fecDG.get());

	GridFunction solution(&fesDG);
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
	ASSERT_NO_THROW(Mesh::LoadFromFile((mfemMeshes2DFolder() + "twotriang.mesh").c_str(), 1, 0));
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

TEST_F(MeshTest, MeshIdentifyBoundaryVertex)
{
	auto mesh{ Mesh::MakeCartesian1D(2) };

	Vector bdrIndex(mesh.GetNBE());
	for (int i = 0; i < bdrIndex.Size(); ++i) {
		bdrIndex(i) = int(mesh.GetBdrElementFaceIndex(i));
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
		bdrIndex(i) = int(mesh.GetBdrElementFaceIndex(i));
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
	auto m{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "FiveQuadsOnX_TFSF.mesh").c_str(), 1, 0) };
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

TEST_F(MeshTest, ExtendedFindNeighboursMethod2D)
{
	auto m{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "square3x3.mesh").c_str(), 1, 0) };
	
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
	auto mesh{ Mesh::LoadFromFile((mfemMeshes1DFolder() + "line.mesh").c_str(), 1, 0) };

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
	auto m{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "four_quads_for_submeshing.mesh").c_str(),1,0) };

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

TEST_F(MeshTest, IdentifyingQuadraticElements_2D)
{
	auto m{ Mesh::LoadFromFile(gmshMeshesFolder() + "2D_Simple_Quadratic.msh",1,0) };
	auto m_copy{ Mesh(m) };
	auto tol{ 1e-3 };

	auto m_fes = m.GetNodalFESpace();
	ASSERT_NE(m_fes->GetNVDofs(), m_fes->GetNDofs());
	for (auto e{ 0 }; e < m.GetNE(); e++) {
		ASSERT_EQ(m.GetElement(0)->GetGeometryType(), m.GetElement(e)->GetGeometryType());
	}

	auto mvol = m.GetElementVolume(0);
	auto mc_vol = m_copy.GetElementVolume(0);
	ASSERT_NEAR(mvol, mc_vol, tol);

	m_copy.SetCurvature(1);

	auto m_copy_fes = m_copy.GetNodalFESpace();
	ASSERT_EQ(m_copy_fes->GetNVDofs(), m_copy_fes->GetNDofs());

	m_copy.SetCurvature(2);

	mc_vol = m_copy.GetElementVolume(0);
	ASSERT_GT(mvol - mc_vol, tol);
	
	ASSERT_EQ(m.GetNodes()->Size(), m_copy.GetNodes()->Size());
	ASSERT_GT(m.GetNodes()->DistanceSquaredTo(*m_copy.GetNodes()), tol);

}

TEST_F(MeshTest, IdentifyingMeshOrder_2D)
{
	auto m_quad { Mesh::LoadFromFile(gmshMeshesFolder() + "2D_Simple_Quadratic.msh", 1, 0) };
	auto m_lin { Mesh::LoadFromFile(gmshMeshesFolder() + "twosquares.msh", 1, 0) };

	EXPECT_EQ(2, m_quad.GetNodalFESpace()->GetMaxElementOrder());

	//Linear meshes do not have a valid Nodes object, it is not defined, thus the call to GetNodalFESpace() 
	//will return a null pointer from which we cannot extract order information.
	//If this call returns said null pointer, we will assume the mesh is linear. If the call returns a valid pointer
	//to the Nodes FES, then we will ask for the MaxElementOrder of the mesh.
	
	auto lin_fes_null{ false };
	if (!m_lin.GetNodalFESpace()) {
		lin_fes_null = true;
	}
	EXPECT_EQ(true, lin_fes_null);

}

TEST_F(MeshTest, IdentifyCurvedElements_2D)
{
	auto gmsh_starting_element{ 15 }; //Gmsh element Id is not unique for each dimensional element. 1D edges count as elements for the id tagging. 15 is the ID for the first 2D element in the mesh.
	Array<ElementId> expectedCurvedElements({ 25, 33, 34, 23, 46, 39, 47, 24, 37, 38, 26, 35, 36 });
	std::sort(expectedCurvedElements.begin(), expectedCurvedElements.end());
	for (auto e{ 0 }; e < expectedCurvedElements.Size(); e++) {
		expectedCurvedElements[e] -= gmsh_starting_element;
	}

	auto mesh_p2{ Mesh::LoadFromFile(gmshMeshesFolder() + "2D_CurvedElementsCircle.msh", 1, 0) };
	DG_FECollection fec_p2(2, mesh_p2.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace fes_p2(&mesh_p2, &fec_p2);
	
	auto lists{ ::initCurvedAndLinearElementsLists(fes_p2, buildDoFPositions(fes_p2)) };
	Array<ElementId> curvedElements;
	for (const auto& [k, v] : lists.second) {
		curvedElements.Append(k);
	}

	EXPECT_EQ(expectedCurvedElements, curvedElements);

}