#pragma once

#include "TestMfemHesthavenFunctions.h"
#include "../TestGlobalFunctions.h"

using namespace mfem;

class MFEMHesthaven3D : public ::testing::Test {
protected:

	void SetUp() override 
	{
		mesh_ = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::HEXAHEDRON);
		fec_ = std::make_unique<DG_FECollection>(1, 3, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
	}

	void set3DFES(
		const int order, 
		const int xElem = 1,
		const int yElem = 1, 
		const int zElem = 1)
	{
		mesh_ = Mesh::MakeCartesian3D(xElem, yElem, zElem, Element::Type::HEXAHEDRON);
		fec_ = std::make_unique<DG_FECollection>(order, 3, BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	}

	Mesh mesh_;
	std::unique_ptr<FiniteElementCollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;

	double tol_ = 1e-6;

};

TEST_F(MFEMHesthaven3D, checkDOperator3DO2)
{
	set3DFES(2);

	auto MISCalcOld = 0.5 * buildInverseMassMatrixEigen(fes_) * buildStiffnessMatrixEigen(fes_);

	std::cout << MISCalcOld << std::endl;

	EXPECT_TRUE(MISCalcOld.isApprox(build3DOneElementDMatrix(), tol_));
}
