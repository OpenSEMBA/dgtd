#include <gtest/gtest.h>
#include <math.h>

#include "TestUtils.h"
#include "math/Geometry.h"
#include "math/Calculus.h"

using namespace maxwell;
using namespace mfem;

class GeometryTest : public ::testing::Test {
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

	EXPECT_TRUE(elementsHaveSameOrientation(&seg1, m.GetFace(m.GetBdrElementFaceIndex(bdr1))));
	EXPECT_FALSE(elementsHaveSameOrientation(&seg2, m.GetFace(m.GetBdrElementFaceIndex(bdr2))));

}

TEST_F(GeometryTest, orientation_from_gmsh_mesh)
{
	auto m{ mfem::Mesh::LoadFromFile((gmshMeshesFolder() + "twosquares.msh").c_str(), 1, 0, true) };

	EXPECT_TRUE(elementsHaveSameOrientation(m.GetFace( m.GetBdrElementFaceIndex(0)), m.GetFace(m.GetBdrElementFaceIndex(1) )));

	EXPECT_FALSE(elementsHaveSameOrientation(m.GetBdrElement(0), m.GetBdrElement(1)));
}

TEST_F(GeometryTest, orientation_from_gmsh_mesh_3D)
{
	//We expect BdrFace and its corresponding face to have opposing normals.

	auto m { Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "test_normals_3D_all_bdr_normals_in.msh").c_str(), 1, 0, false) };

	{ //Bdr Element 2 & Corresponding Face 53 (XY plane, bottom hexa face)
		auto be{ m.GetBdrElement(2) };
		auto be_trans{ m.GetBdrElementTransformation(2) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(2)) };
		EXPECT_EQ(53, m.GetBdrElementFaceIndex(2));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(2)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		for (auto i{ 0 }; i < m.Dimension(); ++i) {
			EXPECT_EQ(normal_be[i], normal_f[i]);
		}
	}

	{ //Bdr Element 4 & Corresponding Face 49 (XZ plane, front hexa face)
		auto be{ m.GetBdrElement(4) };
		auto be_trans{ m.GetBdrElementTransformation(4) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(4)) };
		EXPECT_EQ(49, m.GetBdrElementFaceIndex(4));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(4)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		for (auto i{ 0 }; i < m.Dimension(); ++i) {
			EXPECT_EQ(normal_be[i], normal_f[i]);
		}

	}

	{ //Bdr Element 8 & Corresponding Face 54 (YZ plane, right hexa face)
		auto be{ m.GetBdrElement(8) };
		auto be_trans{ m.GetBdrElementTransformation(8) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(8)) };
		EXPECT_EQ(54, m.GetBdrElementFaceIndex(8));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(8)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		for (auto i{ 0 }; i < m.Dimension(); ++i) {
			EXPECT_EQ(normal_be[i], normal_f[i]);
		}
	}

	{ //Bdr Element 15 & Corresponding Face 30 (XZ plane, back hexa face)
		auto be{ m.GetBdrElement(15) };
		auto be_trans{ m.GetBdrElementTransformation(15) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(15)) };
		EXPECT_EQ(30, m.GetBdrElementFaceIndex(15));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(15)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		for (auto i{ 0 }; i < m.Dimension(); ++i) {
			EXPECT_EQ(normal_be[i], normal_f[i]);
		}
	}

	{ //Bdr Element 19 & Corresponding Face 24 (YZ plane, left hexa face)
		auto be{ m.GetBdrElement(19) };
		auto be_trans{ m.GetBdrElementTransformation(19) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(19)) };
		EXPECT_EQ(24, m.GetBdrElementFaceIndex(19));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(19)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		for (auto i{ 0 }; i < m.Dimension(); ++i) {
			EXPECT_EQ(normal_be[i], normal_f[i]);
		}
	}

	{ //Bdr Element 23 & Corresponding Face 28 (XY plane, top hexa face)
		auto be{ m.GetBdrElement(23) };
		auto be_trans{ m.GetBdrElementTransformation(23) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(23)) };
		EXPECT_EQ(28, m.GetBdrElementFaceIndex(23));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(23)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		for (auto i{ 0 }; i < m.Dimension(); ++i) {
			EXPECT_EQ(normal_be[i], normal_f[i]);
		}
	}

}

TEST_F(GeometryTest, orientation_from_gmsh_mesh_3D_bad_orientation)
{
	//We expect BdrFace and its corresponding face to have opposing normals. 
	//This mesh purposely should fail on surfaces 1, 3, 5, 6 as it lacks a surface reversal. (Gmesh notation) (bottom, right, left, top).

	auto m{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "test_normals_3D_no_reversal.msh").c_str(), 1, 0, false) };

	{ //Bdr Element 2 & Corresponding Face 53 (XY plane, bottom hexa face)
		auto be{ m.GetBdrElement(2) };
		auto be_trans{ m.GetBdrElementTransformation(2) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(2)) };
		EXPECT_EQ(53, m.GetBdrElementFaceIndex(2));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(2)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		EXPECT_FALSE(abs(normal_be.Sum() - normal_f.Sum()) <= 1e-5);
	}

	{ //Bdr Element 4 & Corresponding Face 49 (XZ plane, front hexa face)
		auto be{ m.GetBdrElement(4) };
		auto be_trans{ m.GetBdrElementTransformation(4) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(4)) };
		EXPECT_EQ(49, m.GetBdrElementFaceIndex(4));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(4)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		EXPECT_TRUE(abs(normal_be.Sum() - normal_f.Sum()) <= 1e-5);

	}

	{ //Bdr Element 8 & Corresponding Face 54 (YZ plane, right hexa face)
		auto be{ m.GetBdrElement(8) };
		auto be_trans{ m.GetBdrElementTransformation(8) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(8)) };
		EXPECT_EQ(54, m.GetBdrElementFaceIndex(8));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(8)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		EXPECT_FALSE(abs(normal_be.Sum() - normal_f.Sum()) <= 1e-5);
	}

	{ //Bdr Element 15 & Corresponding Face 30 (XZ plane, back hexa face)
		auto be{ m.GetBdrElement(15) };
		auto be_trans{ m.GetBdrElementTransformation(15) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(15)) };
		EXPECT_EQ(30, m.GetBdrElementFaceIndex(15));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(15)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		EXPECT_TRUE(abs(normal_be.Sum() - normal_f.Sum()) <= 1e-5);
	}

	{ //Bdr Element 19 & Corresponding Face 24 (YZ plane, left hexa face)
		auto be{ m.GetBdrElement(19) };
		auto be_trans{ m.GetBdrElementTransformation(19) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(19)) };
		EXPECT_EQ(24, m.GetBdrElementFaceIndex(19));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(19)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		EXPECT_FALSE(abs(normal_be.Sum() - normal_f.Sum()) <= 1e-5);
	}

	{ //Bdr Element 23 & Corresponding Face 28 (XY plane, top hexa face)
		auto be{ m.GetBdrElement(23) };
		auto be_trans{ m.GetBdrElementTransformation(23) };
		Vector normal_be(3);
		CalcOrtho(be_trans->Jacobian(), normal_be);

		auto face{ m.GetFace(m.GetBdrElementFaceIndex(23)) };
		EXPECT_EQ(28, m.GetBdrElementFaceIndex(23));
		auto f_trans{ m.GetFaceTransformation(m.GetBdrElementFaceIndex(23)) };
		Vector normal_f(3);
		CalcOrtho(f_trans->Jacobian(), normal_f);

		normal_be *= -1.0;
		EXPECT_FALSE(abs(normal_be.Sum() - normal_f.Sum()) <= 1e-5);
	}

}


