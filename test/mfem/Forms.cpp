#include <gtest/gtest.h>

#include <mfem.hpp>
#include <mfemExtension/BilinearIntegrators.h>
#include <mfemExtension/LinearIntegrators.h>
#include <TestUtils.h>
#include <components/RCSManager.cpp>

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

TEST_F(FormTest, LinearForm_w_High_VDIM)
{
	auto vdim = 6;
	auto m = Mesh::MakeCartesian3D(1, 1, 1, Element::Type::HEXAHEDRON);
	auto fec = DG_FECollection(1, m.Dimension());
	auto fes = FiniteElementSpace(&m, &fec, vdim, Ordering::byNODES);

	ConstantCoefficient one(1.0);
	Vector vecone(vdim);
	VectorConstantCoefficient vone(vecone);

	//auto bf{ BilinearForm(&fes) };
	//bf.AddBdrFaceIntegrator(new VectorDivergenceIntegrator(one));
	//bf.Assemble();
	//bf.Finalize();

	auto lf{ LinearForm(&fes) };
	lf = 0.0;
	lf.AddBdrFaceIntegrator(new DGElasticityDirichletLFIntegrator(vone, one, one, 1.0, 1.0));
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
	lf.AddBdrFaceIntegrator(new FarFieldBdrFaceIntegrator(fc, X), bdr_marker);
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
	lf.AddBdrFaceIntegrator(new FarFieldBdrFaceIntegrator(fc, X), bdr_marker);
	lf.Assemble();

}