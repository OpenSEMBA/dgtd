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
	auto mesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };

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
	auto mesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };

	auto bary{ getBarycenterOfFaceElement(mesh, mesh.GetBdrElementFaceIndex(0)) };
	Vector expectedBarycenter({ 1.0, 0.5 });

	for (auto i{ 0 }; i < bary.Size(); i++) {
		EXPECT_NEAR(expectedBarycenter[i], bary[i], tol_);
	}
}

TEST_F(SubMesherTest, barycentersOfTwoElements_2D)
{
	auto mesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };

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
	auto mesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };

	auto baryVector{ buildBarycenterPosition(mesh, 0) };
	Vector expectedBaryVector({ (1.0 + 1.0 / 3.0) - (2.0 / 3.0), (1.0 / 3.0) - (2.0 / 3.0) });

	for (auto i{ 0 }; i < baryVector.Size(); i++) {
		EXPECT_NEAR(expectedBaryVector[i], baryVector[i], tol_);
	}
}

TEST_F(SubMesherTest, tangentVector_2D)
{
	auto mesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };

	auto tangentVector{ buildTangent2D(mesh, 0) };
	Vector expectedTangent({ 0.0, 1.0 });

	for (auto i{ 0 }; i < tangentVector.Size(); i++) {
		EXPECT_NEAR(expectedTangent[i], tangentVector[i], tol_);
	}
}

TEST_F(SubMesherTest, crossProductBaryTangent_2D)
{
	auto mesh{ Mesh::LoadFromFileNoBdrFix(gmshMeshesFolder() + "2D_TwoTriangles_InteriorBdr.msh",1,0) };

	auto faceTrans{ mesh.GetInternalBdrFaceTransformations(0) };
	auto crossValue{ calculateCrossBaryVertexSign(mesh, *faceTrans , 0) };
	auto expectedCrossValue{ 2.0 / 3.0 };

	EXPECT_NEAR(expectedCrossValue, crossValue, tol_);
}

TEST_F(SubMesherTest, calculateNormal_3D)
{
	auto mesh{ Mesh::MakeCartesian3D(1,1,1,Element::Type::HEXAHEDRON) };

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