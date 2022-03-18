#include "FE_Evolution.h"

namespace Maxwell1D {

FE_Evolution::FE_Evolution(FiniteElementSpace* fes) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	fes_(fes),
	MInv_(buildInverseMassMatrix()),
	SE_(buildDerivativeOperator(X, Electric)),
	SH_(buildDerivativeOperator(X, Magnetic)),
	FE_(buildFluxOperator(X, Electric)),
	FH_(buildFluxOperator(X, Magnetic))
{}

std::unique_ptr<BilinearForm> FE_Evolution::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}

std::unique_ptr<BilinearForm> FE_Evolution::buildDerivativeOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto S = std::make_unique<BilinearForm>(fes_);

	ConstantCoefficient one(1.0);

	S->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, d)));
	S->Assemble();
	S->Finalize();

	return S;
}

std::unique_ptr<BilinearForm> FE_Evolution::buildFluxOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto F = std::make_unique<BilinearForm>(fes_);

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	double alpha;
	double beta;

	if (ft == Electric)
	{
		alpha = -1.0;
		beta = 0.0;
		F->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		F->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}
	else
	{
		alpha = -1.0;
		beta = 0.0;
		F->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		F->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}

	F->Assemble();
	F->Finalize();

	return F;
}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(), fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	Vector auxE(MInv_->Height());
	Vector auxFlux(MInv_->Height());
	Vector auxRHS(MInv_->Height());

	FE_->Mult(eOld, auxE);
	FH_->Mult(hOld, auxFlux);
	auxFlux.Add(1.0, auxE);

	// Update E. dE/dt = M^{-1} * (-S * H + ([H] + [E])).
	SH_->AddMult(hOld, auxRHS, -1.0);
	auxRHS.Add(1.0, auxFlux);
	MInv_->Mult(auxRHS, eNew);

	// Update H. dH/dt = M^{-1} * (-S * E + ([E] + [H])).
	SE_->AddMult(eOld, auxRHS, -1.0);
	auxRHS.Add(1.0, auxFlux);
	MInv_->Mult(auxRHS, hNew);
}

}

