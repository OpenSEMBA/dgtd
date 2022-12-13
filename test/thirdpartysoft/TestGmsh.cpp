#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>

using namespace mfem;

class TestGmsh : public ::testing::Test {
	
};

TEST_F(TestGmsh, meshDataGmshMSHRead)
{
	/* This mesh includes a physical tag 1 a surface based on a square split in four trianges with a
	point in the center and a physical tag 2 for the four boundary segments.*/
	ASSERT_NO_THROW(Mesh::LoadFromFile("./TestData/gmshphysicaltags.msh", 1, 0));
}

TEST_F(TestGmsh, 2DboxwithGmshMesh)
{
	auto mesh = Mesh::LoadFromFile("./TestData/gmshphysicaltags.msh", 1, 0);
	auto fec = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
	auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get(), 1, 0);
}

