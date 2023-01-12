#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>
#include "GlobalFunctions.h"
#include "maxwell/mfemExtension/BilinearIntegrators.h"

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

class BilinearIntegratorsExtensionTest : public ::testing::Test {
protected:

	void SetUp() override
	{
		mesh_ = Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE);
		fec_ = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get(), 1, 0);
	}

	void setFES1D(const int order, const int elements = 1)
	{
		mesh_ = Mesh::MakeCartesian1D(elements);
		fec_ = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void setFES2D(const int order, const int nx = 1, const int ny = 1, const double sx = 1.0, const double sy = 1.0, const int vdim = 1)
	{
		mesh_ = Mesh::MakeCartesian2D(nx, ny, Element::Type::TRIANGLE, false, sx, sy);
		fec_ = std::make_unique<DG_FECollection>(order, 2, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get(), vdim, 0);
	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};

TEST_F(BilinearIntegratorsExtensionTest, DISABLED_compareDGTraceWithMaxwellDG1D)
{
	setFES1D(1, 3);

	BilinearForm DGmat(fes_.get());
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], 0.0, 1.0));
	DGmat.Assemble();
	DGmat.Finalize();

	BilinearForm maxwellMat(fes_.get());
	maxwellMat.AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator({X}, 1.0));
	maxwellMat.Assemble();
	maxwellMat.Finalize();

	std::cout << toEigen(*DGmat.SpMat().ToDenseMatrix()) << std::endl;
	std::cout << toEigen(*maxwellMat.SpMat().ToDenseMatrix()) << std::endl;

	EXPECT_TRUE(false);
}


