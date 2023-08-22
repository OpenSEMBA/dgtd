#include <gtest/gtest.h>
#include <math.h>

#include "math/Geometry.h"

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

	EXPECT_TRUE(elementsHaveSameOrientation(&seg1, m.GetFace(m.GetBdrFace(bdr1))));
	EXPECT_FALSE(elementsHaveSameOrientation(&seg2, m.GetFace(m.GetBdrFace(bdr2))));

}