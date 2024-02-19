#include <gtest/gtest.h>
#include <mfem.hpp>
#include <mfemExtension/BilinearIntegrators.h>
#include <mfemExtension/LinearIntegrators.h>
#include <TestUtils.h>
#include <math.h>
#include <complex>

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

double func_exp_real_part_2D(const Vector& x, const double freq, const double phi)
{
	auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
	return cos(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
}

double func_exp_imag_part_2D(const Vector& x, const double freq, const double phi)
{
	auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
	return sin(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
}

class FormTest : public ::testing::Test 
{
public:

	FunctionCoefficient buildFC_2D(const double freq, const double& phi, bool isReal)
	{
		std::function<double(const Vector&)> f = 0;
		switch (isReal) {
		case true:
			f = std::bind(&func_exp_real_part_2D, std::placeholders::_1, freq, phi);
			break;
		case false:
			f = std::bind(&func_exp_imag_part_2D, std::placeholders::_1, freq, phi);
			break;
		}
		FunctionCoefficient res(f);
		return res;
	}
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

TEST_F(FormTest, RCSForm_2D)
{
	
	auto m{ Mesh::LoadFromFile((gmshMeshesFolder() + "2D_Square_RCSTEST.msh").c_str(), 1, 0)};
	auto fec{ DG_FECollection(1, 2) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	GridFunction gf(&fes);
	gf = 2.0;
	LinearForm lf(&fes);
	FunctionCoefficient fc(buildFC_2D(1e6, 0.0, true));
	lf.AddBdrFaceIntegrator(new RCSBdrFaceIntegrator(fc, X));
	lf.Assemble();

}