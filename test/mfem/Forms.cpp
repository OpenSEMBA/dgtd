#include <gtest/gtest.h>

#include <mfem.hpp>
#include <mfemExtension/BilinearIntegrators.h>
#include <mfemExtension/LinearIntegrators.h>
#include <TestUtils.h>
#include <math.h>
#include <complex>
#include <components/RCSManager.cpp>
#include <ctime>

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

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
	
	auto m{ Mesh::LoadFromFile((gmshMeshesFolder() + "2D_LinearForm_RCS.msh").c_str(), 1, 0)};
	auto fec{ DG_FECollection(3, 2) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	GridFunction gf(&fes);
	gf = 2.0;
	LinearForm lf(&fes);
	FunctionCoefficient fc(buildFC_2D(1e7, 0.0, true));
	Array<int> bdr_marker(3);
	bdr_marker = 0;
	bdr_marker[1] = 1;
	lf.AddBdrFaceIntegrator(new RCSBdrFaceIntegrator(fc, X), bdr_marker);
	lf.Assemble();

	auto solution = mfem::InnerProduct(lf, gf);

}

TEST_F(FormTest, RCSBdrFaceInt)
{
	auto m{ Mesh::LoadFromFile((gmshMeshesFolder() + "2D_BdrIntegratorHalfSize.msh").c_str(), 1, 0) };
	auto fec{ DG_FECollection(2,2,BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	LinearForm lf(&fes);
	FunctionCoefficient fc(buildFC_2D(1e7, 0.0, true));
	Array<int> bdr_marker(3);
	bdr_marker = 0;
	bdr_marker[bdr_marker.Size() - 1] = 1;
	lf.AddBdrFaceIntegrator(new RCSBdrFaceIntegrator(fc, X), bdr_marker);
	lf.Assemble();

}

TEST_F(FormTest, disabled_MassMatrixToSingleBlockMatrix)
{
	auto dim = 2;
	auto order = 1;
	auto element_type = Element::TRIANGLE;
	auto nx = 2;
	auto ny = 1;
	auto ne = 0;
	auto dof_per_order = 0;

	switch (element_type) {
	case (Element::Type::TRIANGLE):
		ne = 2 * nx * ny;
		dof_per_order = (order + 1) * (order + 2) / 2;
		break;
	case (Element::Type::QUADRILATERAL):
		ne = nx * ny;
		dof_per_order = int(std::pow(order + 1, 2));
		break;
	}

	Vector numbers(ne * dof_per_order);
	for (int n{ 0 }; n < ne * dof_per_order; n++) {
		numbers[n] = n;
	}

	Vector old_method_result(ne * dof_per_order);

	clock_t full_matrix_clock_start{ clock() };

	auto m{ Mesh::MakeCartesian2D(nx, ny, element_type, false, 2.0, 2.0) };
	auto fec{ DG_FECollection(order, dim, BasisType::GaussLegendre) };
	auto fes{ FiniteElementSpace(&m, &fec) };

	BilinearForm mass(&fes);
	ConstantCoefficient coeff_one(1.0);
	mass.AddDomainIntegrator(new MassIntegrator(coeff_one));
	mass.Assemble();
	mass.Finalize();

	for (auto rep{ 0 }; rep < 100000; rep++) {
		mass.Mult(numbers, old_method_result);
	}

	clock_t full_matrix_clock_end{ clock() };

	double full_matrix_clock_elapsed{ static_cast<double>(full_matrix_clock_end - full_matrix_clock_start) / CLOCKS_PER_SEC };

	for (auto e{ 0 }; e < m.GetNE(); e++) {
		assert(element_type == m.GetElementType(e));
	}

	clock_t single_matrix_clock_start{ clock() };

	auto m_single{ Mesh::MakeCartesian2D(1, 1, element_type) };
	auto fes_single{ FiniteElementSpace(&m_single, &fec) };

	BilinearForm mass_single(&fes_single);
	mass_single.AddDomainIntegrator(new MassIntegrator(coeff_one));
	mass_single.Assemble();
	mass_single.Finalize();

	Vector reduced_mass_coefficients(ne);
	reduced_mass_coefficients = 0.0;
	for (auto d{ 0 }; d < ne; d++) {
		reduced_mass_coefficients[d] = mass.SpMat().ToDenseMatrix()->Elem(d + dof_per_order, d + dof_per_order) /
			mass_single.SpMat().ToDenseMatrix()->Elem(0, 0);
	}

	Vector new_method_result(ne * dof_per_order);
	auto dense = mass_single.SpMat().ToDenseMatrix();
	for (auto rep{ 0 }; rep < 100000; rep++) {
	new_method_result = 0.0;
	for (auto e{ 0 }; e < ne; e++) {
		for (auto d{ 0 }; d < dof_per_order; d++) {
			for (auto d2{ 0 }; d2 < dof_per_order; d2++) {
				new_method_result[e * dof_per_order + d] +=
					reduced_mass_coefficients[e] *
					dense->Elem(d, d2) *
					numbers[e * dof_per_order + d2];
				}
			}
		}
	}

	clock_t single_matrix_clock_end{ clock() };

	double single_matrix_clock_elapsed{ static_cast<double>(single_matrix_clock_end - single_matrix_clock_start) / CLOCKS_PER_SEC };

	std::cout << "Program for full matrix in " << full_matrix_clock_elapsed << " seconds." << std::endl;
	std::cout << "Program for single matrix in " << single_matrix_clock_elapsed << " seconds." << std::endl;
	
	double tol{ 1e-2 };
	for (int v{ 0 }; v < new_method_result.Size(); v++) {
		ASSERT_NEAR(new_method_result[v], old_method_result[v], tol);
	}

	ASSERT_TRUE(false);

}