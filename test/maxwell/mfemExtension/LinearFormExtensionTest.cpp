#include <gtest/gtest.h>

#include "mfemExtension/LinearIntegrators.h"
#include "components/Types.h"

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

class LinearFormExtensionTest : public ::testing::Test 
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

void TDgaussian_function(const Vector& pos, double time, Vector& v)
{
	switch (pos.Size()) {
	case 1:
		v(0) = 1.0 *
			exp(
				-pow(pos[X] - time, 2) /
				(2.0 * pow(0.3, 2))
			);
		break;
	default:
		throw std::runtime_error("Invalid dimension.");
	}
};

class SimpleFEEvol : public TimeDependentOperator
{
private:
	BilinearForm& M_;
	LinearForm& b_;

public:
	SimpleFEEvol(BilinearForm& M, LinearForm& b);

	virtual void Mult(const Vector& x, Vector& y) const;
};

SimpleFEEvol::SimpleFEEvol(BilinearForm& M, LinearForm& b) : TimeDependentOperator(M.FESpace()->GetNDofs()),
M_{ M },
b_{ b }
{}

void SimpleFEEvol::Mult(const Vector& x, Vector& y) const
{
	//b_.Assemble();
	Vector temp(x.Size());
	add(x, b_, temp);
	M_.Mult(temp, y);
}
TEST_F(LinearFormExtensionTest, checkLinearFormFunctionUsage)
{
	int order{ 2 }, dim{ 1 };
	Mesh mesh{ Mesh::MakeCartesian1D(3,1.0)};
	DG_FECollection fec{ order, dim, BasisType::GaussLobatto };
	FiniteElementSpace fes{ &mesh, &fec,1,Ordering::byNODES };

	VectorConstantCoefficient one(Vector({ 1.0 }));

	BilinearForm m{&fes};
	m.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));

	LinearForm b{ &fes };
	b.AddBdrFaceIntegrator(
		new BoundaryDGJumpIntegrator{ one, 1.0 }
	);
	
	m.Assemble();
	m.Finalize();
	b.Assemble();
	
	std::unique_ptr<ODESolver> odeSolver = std::make_unique<RK4Solver>() ;

	std::unique_ptr<TimeDependentOperator> evol = std::make_unique<SimpleFEEvol>(m, b);
	
	double time = 0.0;
	double tf = 1.0;
	double timeStep = 5e-3;

	evol.get()->SetTime(time);
	odeSolver->Init(*evol);

	GridFunction field(&fes);
	field = 0.0;

	for (time; time <= tf; time += timeStep)
	{
		//tdGaussian.SetTime(time);
		//evol.get()->SetTime(time);
		double dt_real = std::min(timeStep, tf - time);
		odeSolver->Step(field, time, dt_real);
	}

}

