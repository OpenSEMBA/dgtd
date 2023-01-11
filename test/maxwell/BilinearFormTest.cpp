#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>
#include "GlobalFunctions.h"
#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/mfemExtension/BilinearForm.h"

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

class BilinearFormTest : public ::testing::Test 
{
};

TEST_F(BilinearFormTest, bdrInteriorBilinearForm)
{
	auto mesh{ Mesh::MakeCartesian1D(4,4.0) };
	auto bdrIntAtt{ 5 };
	mesh.AddBdrPoint(2, bdrIntAtt);
	mesh.FinalizeMesh();

	auto fec{ DG_FECollection(1,1,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&mesh,&fec) };

	BilinearFormTF totalFieldFlux{ &fes };
	Array<int> bdrMarker{ mesh.bdr_attributes.Max() };
	bdrMarker = 0;
	bdrMarker[bdrIntAtt - 1] = 1;
	totalFieldFlux.AddBdrFaceIntegrator(
		new mfemExtension::MaxwellDGTraceJumpIntegrator{
			std::vector<Direction>{X}, 1.0
		},
		bdrMarker
	);
	totalFieldFlux.Assemble();
	totalFieldFlux.Finalize();

	GridFunction f{ &fes }, exc{ &fes };
	f = 0.0; exc = 3.0;
	
	totalFieldFlux.Mult(exc, f);

	EXPECT_EQ(0.0, f[0]);
	EXPECT_EQ(0.0, f[7]);
	EXPECT_EQ(3.0, f[3]);
}