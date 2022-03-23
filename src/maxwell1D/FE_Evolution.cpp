#include "FE_Evolution.h"

namespace Maxwell1D {

FE_Evolution::FE_Evolution(FiniteElementSpace* fes) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	fes_(fes),
	MInv_(buildInverseMassMatrix())
{
	KEE_ = std::make_unique<BilinearForm>(fes_);
	KEH_ = std::make_unique<BilinearForm>(fes_);
	KHE_ = std::make_unique<BilinearForm>(fes_);
	KHH_ = std::make_unique<BilinearForm>(fes_);
	initializeBilinearForms();
}

std::unique_ptr<BilinearForm> FE_Evolution::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}

void FE_Evolution::addDerivativeOperator(
	std::unique_ptr<BilinearForm>& form,
	const Direction& d, ConstantCoefficient& coeff) const
{
	assert(d == X, "Incorrect argument for direction.");

	switch (fluxType) {
		case Upwind:
			form->AddDomainIntegrator(
				new DerivativeIntegrator(coeff, d));
			break;
		case Centered:
			form->AddDomainIntegrator(
				new TransposeIntegrator(
					new DerivativeIntegrator(coeff, d)));
			break;
	}
}

void FE_Evolution::addFluxOperator(
	std::unique_ptr<BilinearForm>& form,
	const Direction& d, Vector& abg) const
{
	assert(d == X, "Incorrect argument for direction.");

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	switch (fluxType) {
		case Upwind:
			form->AddInteriorFaceIntegrator(
				new MaxwellDGTraceIntegrator(n[d], abg[Alpha], abg[Beta], abg[Gamma]));
			form->AddBdrFaceIntegrator(
				new MaxwellDGTraceIntegrator(n[d], abg[Alpha], abg[Beta], abg[Gamma]));
			break;

		case Centered:
			form->AddInteriorFaceIntegrator(
				new DGTraceIntegrator(n[d], abg[Alpha], abg[Beta]));
			form->AddBdrFaceIntegrator(
				new DGTraceIntegrator(n[d], abg[Alpha], abg[Beta]));
			break;

	}
}


void FE_Evolution::finalizeBilinearForm(std::unique_ptr<BilinearForm>& form) const
{
	form->Assemble();
	form->Finalize();
}

void FE_Evolution::initializeBilinearForms()
{
	double alpha = 0.0;
	double beta = 0.0;
	double gamma = 0.0;
	ConstantCoefficient coeff;
	Vector abg({alpha, beta, gamma});

	switch (fluxType) {
	case Upwind:

		abg[Gamma] = 1.0;

		addFluxOperator(KEE_, X, abg);
		finalizeBilinearForm(KEE_);

		addFluxOperator(KHH_, X, abg);
		finalizeBilinearForm(KHH_);

		abg[Gamma] = 1.0 * 0.5;
		coeff = ConstantCoefficient(-1.0);

		addDerivativeOperator(KEH_, X, coeff);
		addFluxOperator(KEH_, X, abg);
		finalizeBilinearForm(KEH_);

		addDerivativeOperator(KHE_, X, coeff);
		addFluxOperator(KHE_, X, abg);
		finalizeBilinearForm(KHE_);

	case Centered:

		abg[Alpha] = -1.0;

		coeff = ConstantCoefficient(-1.0);

		addDerivativeOperator(KEH_, X, coeff);
		addFluxOperator(KEH_, X, abg);
		finalizeBilinearForm(KEH_);

		addDerivativeOperator(KHE_, X, coeff);
		addFluxOperator(KHE_, X, abg);
		finalizeBilinearForm(KHE_);

	}
}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(), fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	Vector auxFluxdE(MInv_->Height());
	Vector auxFluxdH(MInv_->Height());
	Vector auxRHS(MInv_->Height());

	switch (fluxType) {
		case Upwind:

			// Update E. dE/dt = M^{-1} * (-S * H + 0.5 * ([H] - [E])).
			KEH_->Mult(hOld, auxRHS);
			KEE_->Mult(eOld, auxFluxdE);
			auxRHS.Add(-1.0, auxFluxdE);
			MInv_->Mult(auxRHS, eNew);

			// Update H. dH/dt = M^{-1} * (-S * E + 0.5 * ([E] - [H])).
			KHE_->Mult(eOld, auxRHS);
			KHH_->Mult(hOld, auxFluxdH);
			auxRHS.Add(-1.0, auxFluxdH);
			MInv_->Mult(auxRHS, eNew);

			break;

		case Centered:
	
			// Update E. dE/dt = M^{-1} * (-S * H + {H}).
			KEH_->Mult(hOld, auxRHS);
			MInv_->Mult(auxRHS, eNew);

			// Update H. dH/dt = M^{-1} * (-S * E + {E}).
			KHE_->Mult(eOld, auxRHS);
			MInv_->Mult(auxRHS, hNew);

			break;
		}

}

}

