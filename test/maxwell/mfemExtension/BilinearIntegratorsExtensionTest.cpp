#include <gtest/gtest.h>

#include "math/EigenMfemTools.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "components/Types.h"
#include "TestUtils.h"

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

TEST_F(BilinearIntegratorsExtensionTest, checkNormalsOnRotatedQuad)
{
	int dim{ 2 };
	auto m{ Mesh::LoadFromFileNoBdrFix((mfemMeshes2DFolder() + "rotatedquad.mesh"),1,0) };
	DG_FECollection fec{ 1, dim, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &m, &fec };

	BilinearForm bf(&fes);
	bf.AddBdrFaceIntegrator(new MaxwellDGTraceJumpIntegrator({ X,Y }, 1.0));
	bf.Assemble();
	bf.Finalize();

}