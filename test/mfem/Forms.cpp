#include <gtest/gtest.h>
#include <mfem.hpp>
#include <mfemExtension/BilinearIntegrators.h>
#include <mfemExtension/LinearIntegrators.h>
#include <TestUtils.h>

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

class FormTest : public ::testing::Test {
};

TEST_F(FormTest, LinearForms_1D) 
{

	auto m{ Mesh::MakeCartesian1D(3) };
	auto att = 301;
	m.AddBdrPoint(1, att);
	m.FinalizeMesh();

	auto fec{ DG_FECollection(1,1,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	auto vec{ Vector(1) }; vec[0] = 5.0;
	auto vcc{ VectorConstantCoefficient(vec) };
	Array<int> int_att_bdr_marker(m.bdr_attributes.Max());
	int_att_bdr_marker[att - 1] = 1;

	auto lf{ LinearForm(&fes) };
	lf = 0.0;
	lf.AddInternalBoundaryFaceIntegrator(
		new BoundaryDGJumpIntegrator(vcc, 1.0), int_att_bdr_marker
	);
	
	lf.Assemble();

}

TEST_F(FormTest, LinearForms_2D) 
{
	
	auto m{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "intbdr_two_quads.mesh").c_str(), 1, 0) };
	auto fec{ DG_FECollection(1,2,BasisType::GaussLegendre) };
	auto fes{ FiniteElementSpace(&m,&fec) };

	auto vec{ Vector(2) }; vec[0] = 5.0; vec[1] = 2.0;
	auto vcc{ VectorConstantCoefficient(vec) };
	Array<int> int_att_bdr_marker(m.bdr_attributes.Max());
	int_att_bdr_marker[m.bdr_attributes.Max() - 1] = 1;

	auto lf{ LinearForm(&fes) };
	lf = 0.0;
	lf.AddInternalBoundaryFaceIntegrator(
		new BoundaryDGJumpIntegrator(vcc, 1.0), int_att_bdr_marker
	);

	lf.Assemble();
	
}