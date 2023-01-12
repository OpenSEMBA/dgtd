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
		const double length = 1.0)
	{
		mesh_ = Mesh::MakeCartesian1D(elements, length);
		fec_ = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

};

TEST_F(BilinearFormTest, bdrInteriorBilinearForm)
{
	setFES1D(1,4,4.0);

	auto intBdrAttr{ 5 };
	mesh_.AddBdrPoint(2, intBdrAttr);
	mesh_.FinalizeMesh();

	BilinearFormTF totalFieldFlux{ fes_.get() };
	Array<int> intBdrMarker{ mesh_.bdr_attributes.Max() };
	intBdrMarker = 0;
	intBdrMarker[intBdrAttr - 1] = 1;
	std::vector<VectorConstantCoefficient> n = {VectorConstantCoefficient(Vector({1.0}))};
	totalFieldFlux.AddInteriorBoundaryFaceIntegrator(
		new DGTraceIntegrator{n[0],0.0, 1.0},
		intBdrMarker
	);
	totalFieldFlux.Assemble();
	totalFieldFlux.Finalize();

	GridFunction f{ fes_.get() }, exc{ fes_.get() };
	f = 0.0; 
	exc[3] = 6.0;
	exc[4] = 5.0;
	
	totalFieldFlux.Mult(exc, f);

	EXPECT_EQ( 0.0, f[0]);
	EXPECT_EQ( 0.0, f[7]);
	EXPECT_EQ( 1.0, f[3]);
	EXPECT_EQ(-1.0, f[4]);
}