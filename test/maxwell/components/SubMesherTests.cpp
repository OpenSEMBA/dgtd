#include "components/SubMesher.h"
#include "TestUtils.h"

#include <gtest/gtest.h>
#include <math.h>

using namespace maxwell;
using namespace mfem;

class SubMesherTest : public ::testing::Test {
public:
	double tol_ = 1e-4;
};

TEST_F(SubMesherTest, barycenterOfElements_2D)
{
	auto smesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	{
		auto bary{ getBarycenterOfElement(mesh, 0) };
		Vector expectedBarycenter({ 1.0 / 3.0, 1.0 / 3.0 });
		for (auto i{ 0 }; i < bary.Size(); i++) {
			EXPECT_NEAR(expectedBarycenter[i], bary[i], tol_);
		}
	}

	{
		auto bary{ getBarycenterOfElement(mesh, 1) };
		Vector expectedBarycenter({ 2.0 / 3.0, 2.0 / 3.0 });
		for (auto i{ 0 }; i < bary.Size(); i++) {
			EXPECT_NEAR(expectedBarycenter[i], bary[i], tol_);
		}
	}

	{
		auto bary{ getBarycenterOfElement(mesh, 2) };
		Vector expectedBarycenter({ 1.0 + 1.0 / 3.0, 1.0 / 3.0 });
		for (auto i{ 0 }; i < bary.Size(); i++) {
			EXPECT_NEAR(expectedBarycenter[i], bary[i], tol_);
		}
	}

	{
		auto bary{ getBarycenterOfElement(mesh, 3) };
		Vector expectedBarycenter({ 1.0 + 2.0 / 3.0, 2.0 / 3.0 });
		for (auto i{ 0 }; i < bary.Size(); i++) {
			EXPECT_NEAR(expectedBarycenter[i], bary[i], tol_);
		}
	}
}

TEST_F(SubMesherTest, barycenterOfFaceElements_2D)
{
	auto smesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	auto bary{ getBarycenterOfFaceElement(mesh, mesh.GetBdrElementFaceIndex(0)) };
	Vector expectedBarycenter({ 1.0, 0.5 });

	for (auto i{ 0 }; i < bary.Size(); i++) {
		EXPECT_NEAR(expectedBarycenter[i], bary[i], tol_);
	}
}

TEST_F(SubMesherTest, barycentersOfTwoElements_2D)
{
	auto smesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	auto barycenters{ calculateBarycenters(mesh, 0) };
	Vector expectedFirstBarycenter({ 2.0 / 3.0, 2.0 / 3.0 });
	Vector expectedSecondBarycenter({ 1.0 + 1.0 / 3.0, 1.0 / 3.0 });

	for (auto i{ 0 }; i < expectedFirstBarycenter.Size(); i++) {
		EXPECT_NEAR(expectedFirstBarycenter[i], barycenters.first[i], tol_);
		EXPECT_NEAR(expectedSecondBarycenter[i], barycenters.second[i], tol_);
	}
}

TEST_F(SubMesherTest, barycenterVectorOfTwoElements_2D)
{
	auto smesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	auto baryVector{ buildElem1ToElem2BarycenterVector(mesh, 0) };
	Vector expectedBaryVector({ (1.0 + 1.0 / 3.0) - (2.0 / 3.0), (1.0 / 3.0) - (2.0 / 3.0) });

	for (auto i{ 0 }; i < baryVector.Size(); i++) {
		EXPECT_NEAR(expectedBaryVector[i], baryVector[i], tol_);
	}
}

TEST_F(SubMesherTest, tangentVector_2D)
{
	auto smesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	auto tangentVector{ buildTangent2D(mesh, 0) };
	Vector expectedTangent({ 0.0, 1.0 });

	for (auto i{ 0 }; i < tangentVector.Size(); i++) {
		EXPECT_NEAR(expectedTangent[i], tangentVector[i], tol_);
	}
}

TEST_F(SubMesherTest, crossProductBaryTangent_2D)
{
	auto smesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	auto faceTrans{ mesh.GetInternalBdrFaceTransformations(0) };
	auto crossValue{ calculateCrossBaryVertexSign(mesh, *faceTrans , 0) };
	auto expectedCrossValue{ 2.0 / 3.0 };

	EXPECT_NEAR(expectedCrossValue, crossValue, tol_);
}

TEST_F(SubMesherTest, calculateNormal_3D)
{
	auto smesh{ Mesh::MakeCartesian3D(1,1,1,Element::Type::HEXAHEDRON) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);

	{
		auto normal = buildNormal3D(mesh, 0);
		Vector expectedNormal({ 0.0, 0.0, -1.0 });
		for (auto i{ 0 }; i < normal.Size(); i++) {
			EXPECT_NEAR(expectedNormal[i], normal[i], tol_);
		}
	}
	{
		auto normal = buildNormal3D(mesh, 1);
		Vector expectedNormal({ 0.0, 0.0,  1.0 });
		for (auto i{ 0 }; i < normal.Size(); i++) {
			EXPECT_NEAR(expectedNormal[i], normal[i], tol_);
		}
	}
	{
		auto normal = buildNormal3D(mesh, 2);
		Vector expectedNormal({ -1.0, 0.0, 0.0 });
		for (auto i{ 0 }; i < normal.Size(); i++) {
			EXPECT_NEAR(expectedNormal[i], normal[i], tol_);
		}
	}
	{
		auto normal = buildNormal3D(mesh, 3);
		Vector expectedNormal({  1.0, 0.0, 0.0 });
		for (auto i{ 0 }; i < normal.Size(); i++) {
			EXPECT_NEAR(expectedNormal[i], normal[i], tol_);
		}
	}
	{
		auto normal = buildNormal3D(mesh, 4);
		Vector expectedNormal({ 0.0, -1.0, 0.0 });
		for (auto i{ 0 }; i < normal.Size(); i++) {
			EXPECT_NEAR(expectedNormal[i], normal[i], tol_);
		}
	}
	{
		auto normal = buildNormal3D(mesh, 5);
		Vector expectedNormal({ 0.0,  1.0, 0.0 });
		for (auto i{ 0 }; i < normal.Size(); i++) {
			EXPECT_NEAR(expectedNormal[i], normal[i], tol_);
		}
	}
}

TEST_F(SubMesherTest, TFSFSubMesher1D_OnePointTwoPointsAndThreePointRejection)
{
	// Build a 4-element 1D mesh on [0, 4] with interior boundary points at vertices 1, 2, 3
	Mesh mesh(1, 5, 4, 5);
	for (int i = 0; i < 5; ++i) {
		const double x = static_cast<double>(i);
		mesh.AddVertex(&x);
	}
	for (int e = 0; e < 4; ++e) {
		mesh.AddSegment(e, e + 1, 1);
	}
	auto addBdrPt = [&](int v, int attr) {
		Element* b = mesh.NewElement(Geometry::POINT);
		b->SetVertices(&v);
		b->SetAttribute(attr);
		mesh.AddBdrElement(b);
	};
	addBdrPt(0, 1); // left outer boundary
	addBdrPt(4, 2); // right outer boundary
	addBdrPt(1, 3); // interior point x = 1 (between el 0 and el 1)
	addBdrPt(3, 4); // interior point x = 3 (between el 2 and el 3)
	addBdrPt(2, 5); // interior point x = 2 (between el 1 and el 2)
	mesh.FinalizeTopology();
	mesh.Finalize();

	// 1. Single-point 1D TFSF at tag 3 (x = 1): el 0 is SF, el 1 is TF
	{
		Array<int> marker(5);
		marker = 0;
		marker[3 - 1] = 1;
		TotalFieldScatteredFieldSubMesher submesher(mesh, marker);
		ASSERT_NE(submesher.getSFSubMesh(), nullptr);
		ASSERT_NE(submesher.getTFSubMesh(), nullptr);
		const auto& sf_map = submesher.getSFSubMesh()->GetParentElementIDMap();
		const auto& tf_map = submesher.getTFSubMesh()->GetParentElementIDMap();
		ASSERT_EQ(sf_map.Size(), 1);
		ASSERT_EQ(tf_map.Size(), 1);
		EXPECT_EQ(sf_map[0], 0);
		EXPECT_EQ(tf_map[0], 1);
	}

	// 2. Two-point 1D TFSF at tags 3 (x = 1) and 4 (x = 3):
	// outer elements {0, 3} are SF, inner adjacent elements {1, 2} are TF
	{
		Array<int> marker(5);
		marker = 0;
		marker[3 - 1] = 1;
		marker[4 - 1] = 1;
		TotalFieldScatteredFieldSubMesher submesher(mesh, marker);
		ASSERT_NE(submesher.getSFSubMesh(), nullptr);
		ASSERT_NE(submesher.getTFSubMesh(), nullptr);
		const auto& sf_map = submesher.getSFSubMesh()->GetParentElementIDMap();
		const auto& tf_map = submesher.getTFSubMesh()->GetParentElementIDMap();
		ASSERT_EQ(sf_map.Size(), 2);
		ASSERT_EQ(tf_map.Size(), 2);
		EXPECT_EQ(sf_map[0], 0);
		EXPECT_EQ(sf_map[1], 3);
		EXPECT_EQ(tf_map[0], 1);
		EXPECT_EQ(tf_map[1], 2);
	}

	// 3. Three TFSF points in 1D must throw std::runtime_error
	{
		Array<int> marker(5);
		marker = 0;
		marker[3 - 1] = 1;
		marker[4 - 1] = 1;
		marker[5 - 1] = 1;
		EXPECT_THROW(TotalFieldScatteredFieldSubMesher(mesh, marker), std::runtime_error);
	}
}