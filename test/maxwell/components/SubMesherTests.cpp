#include "components/SubMesher.h"
#include "TestUtils.h"

#include <gtest/gtest.h>
#include <math.h>

using namespace maxwell;
using namespace mfem;

static Mesh make1DPMLShellMesh(int n = 5)
{
	auto m = Mesh::MakeCartesian1D(n);
	m.GetElement(0)->SetAttribute(2);
	m.GetElement(n - 1)->SetAttribute(2);
	m.SetAttributes();
	return m;
}

static Array<int> makePMLMarker()
{
	Array<int> marker(2);
	marker = 0;
	marker[1] = 1;
	return marker;
}

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

TEST_F(SubMesherTest, volumetricPMLSubMesher_extractsSerialShellDomain)
{
	auto mesh = make1DPMLShellMesh();
	auto marker = makePMLMarker();

	VolumetricPMLSubMesher submesher(mesh, marker);
	auto* submesh = submesher.getSubMesh();
	ASSERT_NE(nullptr, submesh);
	ASSERT_EQ(2, submesh->GetNE());
	ASSERT_EQ(2, submesher.getParentElementIDMap().Size());
	EXPECT_EQ(0, submesher.getParentElementIDMap()[0]);
	EXPECT_EQ(4, submesher.getParentElementIDMap()[1]);
	ASSERT_EQ(4, submesher.getParentVertexIDMap().Size());
	EXPECT_EQ(0, submesher.getParentVertexIDMap()[0]);
	EXPECT_EQ(1, submesher.getParentVertexIDMap()[1]);
	EXPECT_EQ(4, submesher.getParentVertexIDMap()[2]);
	EXPECT_EQ(5, submesher.getParentVertexIDMap()[3]);
}

TEST_F(SubMesherTest, volumetricPMLSubMesher_extractsParallelShellDomain)
{
	auto serial_mesh = make1DPMLShellMesh();
	ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
	auto marker = makePMLMarker();

	VolumetricPMLSubMesher submesher(mesh, marker);
	auto* submesh = submesher.getParSubMesh();
	ASSERT_NE(nullptr, submesh);

	const int local_count = submesher.getParentElementIDMap().Size();
	int global_count = 0;
	MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	EXPECT_EQ(2, global_count);

	for (int i = 0; i < submesher.getParentElementIDMap().Size(); ++i) {
		EXPECT_EQ(2, mesh.GetAttribute(submesher.getParentElementIDMap()[i]));
	}

	EXPECT_EQ(local_count + 2, submesher.getParentVertexIDMap().Size());
	for (int i = 0; i < submesher.getParentVertexIDMap().Size(); ++i) {
		EXPECT_LE(0, submesher.getParentVertexIDMap()[i]);
	}
}

TEST_F(SubMesherTest, volumetricPMLSubMesher_requiresMarkedDomain)
{
	auto mesh = make1DPMLShellMesh();
	Array<int> marker(2);
	marker = 0;

	EXPECT_THROW(VolumetricPMLSubMesher(mesh, marker), std::runtime_error);
}