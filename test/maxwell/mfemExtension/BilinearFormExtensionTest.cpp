#include <gtest/gtest.h>

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"
#include "components/Types.h"

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

class BilinearFormExtensionTest : public ::testing::Test 
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
		const double length = 1.0,
		const decltype(BasisType::GaussLobatto) basis = BasisType::GaussLobatto
	)
	{
		mesh_ = Mesh::MakeCartesian1D(elements, length);
		fec_ = std::make_unique<DG_FECollection>(order, 1, basis);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;
};

TEST_F(BilinearFormExtensionTest, checkInteriorBoundaryFaceIntegrator)
{
	setFES1D(1,4,4.0);

	auto intBdrAttr{ 5 };
	mesh_.AddBdrPoint(2, intBdrAttr);
	mesh_.FinalizeMesh();

	BilinearFormIBFI totalFieldFlux{ fes_.get() };
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

TEST_F(BilinearFormExtensionTest, compareBaseAndDerivedBilinearForms)
{
	setFES1D(1, 4, 4.0);

	auto intBdrAttr{ 5 };
	mesh_.AddBdrPoint(2, intBdrAttr);
	mesh_.FinalizeMesh();

	Array<int> intBdrMarker{ mesh_.bdr_attributes.Max() };
	intBdrMarker = 0;
	intBdrMarker[intBdrAttr - 1] = 1;
	std::vector<VectorConstantCoefficient> n = { VectorConstantCoefficient(Vector({1.0})) };

	BilinearForm baseBilinearForm(fes_.get());
	BilinearFormIBFI derivedBilinearForm(fes_.get());
	baseBilinearForm.AddBdrFaceIntegrator(
		new DGTraceIntegrator{ n[0],0.0, 1.0 },
		intBdrMarker
	);
	derivedBilinearForm.AddInteriorBoundaryFaceIntegrator(
		new DGTraceIntegrator{ n[0],0.0, 1.0 },
		intBdrMarker
	);
	baseBilinearForm.Assemble();
	baseBilinearForm.Finalize();
	derivedBilinearForm.Assemble();
	derivedBilinearForm.Finalize();	

	GridFunction fbase{ fes_.get() }, fderived{ fes_.get() },exc{ fes_.get() };
	exc[3] = 6.0;
	exc[4] = 5.0;
	baseBilinearForm.Mult(exc, fbase);
	derivedBilinearForm.Mult(exc, fderived);

	EXPECT_EQ( 0.0, fbase[3]);
	EXPECT_EQ( 0.0, fbase[4]);
	EXPECT_EQ( 1.0, fderived[3]);
	EXPECT_EQ(-1.0, fderived[4]);
	
}


