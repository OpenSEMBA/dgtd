#include "FE_Evolution.h"

namespace Maxwell1D {

FE_Evolution::FE_Evolution(FiniteElementSpace* fes) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	fes_(fes),
	MInv_(buildInverseMassMatrix())
{
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
	const Direction& d, const  Vector& abgFace, const Vector& abgBdr) const
{
	assert(d == X, "Incorrect argument for direction.");

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	switch (fluxType) {
	case Upwind:
		form->AddInteriorFaceIntegrator(
			new MaxwellDGTraceIntegrator(n[d], abgFace[Alpha], abgFace[Beta]));
		form->AddBdrFaceIntegrator(
			new MaxwellDGTraceIntegrator(n[d], abgBdr[Alpha], abgBdr[Beta]));
		break;

	case Centered:
		form->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], abgFace[Alpha], abgFace[Beta]));
		form->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], abgBdr[Alpha], abgBdr[Beta]));
		break;

	}
	form->Assemble(0);
	form->Finalize(0);
}

void FE_Evolution::initializeBilinearForms()
{
	double alpha = 0.0;
	double beta = 0.0;
	ConstantCoefficient coeff;
	Vector abIntFaceE({alpha, beta});
	Vector abBdrFaceE({alpha, beta});
	Vector abIntFaceH({alpha, beta});
	Vector abBdrFaceH({alpha, beta});

	KEE_ = std::make_unique<BilinearForm>(fes_);
	KEH_ = std::make_unique<BilinearForm>(fes_);
	KHE_ = std::make_unique<BilinearForm>(fes_);
	KHH_ = std::make_unique<BilinearForm>(fes_);

	switch (fluxType) {
	case Upwind:

		switch (bdrCond) {
		case PEC:
			abIntFaceE[Alpha] = 0.0;
			abIntFaceH[Alpha] = 1.0;
			abBdrFaceE[Alpha] = 0.0;
			abBdrFaceH[Alpha] = 1.0;
			abIntFaceE[Beta]  = 0.5;
			abIntFaceH[Beta]  = 0.0;
			abBdrFaceE[Beta]  = 0.5;
			abBdrFaceH[Beta]  = 0.0;
			break;		
		default:		
			abIntFaceE[Beta] = -1.0;
			abIntFaceH[Beta] = -1.0;
			abBdrFaceE[Beta] = -1.0;
			abBdrFaceH[Beta] = -1.0;
			break;
		}

		coeff = ConstantCoefficient(-1.0);

		addFluxOperator(KEE_, X, abIntFaceE, abBdrFaceE);

		addDerivativeOperator(KEH_, X, coeff);
		addFluxOperator(KEH_, X, abIntFaceE, abBdrFaceE);

		switch (bdrCond) {
		case PEC:
			abIntFaceE[Alpha] = 1.0;
			abIntFaceH[Alpha] = 0.0;
			abBdrFaceE[Alpha] = 1.0;
			abBdrFaceH[Alpha] = 0.0;
			abIntFaceE[Beta]  = 0.0;
			abIntFaceH[Beta]  = 0.5;
			abBdrFaceE[Beta]  = 0.0;
			abBdrFaceH[Beta]  = 0.5;
			break;
		default:
			throw std::exception("Input a valid Boundary Condition.");
			break;
		}

		addFluxOperator(KHH_, X, abIntFaceH, abBdrFaceH);

		addDerivativeOperator(KHE_, X, coeff);
		addFluxOperator(KHE_, X, abIntFaceH, abBdrFaceH);

		break;

	case Centered:

		switch (bdrCond) {
		case PEC:
			abIntFaceE[Alpha] = -1.0;
			abIntFaceH[Alpha] = -1.0;
			abBdrFaceE[Alpha] =  0.0;
			abBdrFaceH[Alpha] = -2.0;
			break;
		case PMC:
			abIntFaceE[Alpha] = -1.0;
			abIntFaceH[Alpha] = -1.0;
			abBdrFaceE[Alpha] = -2.0;
			abBdrFaceH[Alpha] =  0.0;
			break;
		case SMA:
			throw std::exception("SMA not available for Centered Flux.");
		default:
			throw std::exception("Input a valid Boundary Condition.");
			break;
		}

		coeff = ConstantCoefficient(1.0);

		addDerivativeOperator(KEH_, X, coeff);
		addFluxOperator(KEH_, X, abIntFaceE, abBdrFaceE);

		addDerivativeOperator(KHE_, X, coeff);
		addFluxOperator(KHE_, X, abIntFaceH, abBdrFaceH);

		break;
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

		// Update E. dE/dt = M^{-1} * (-S * H + nx * {H} + 0.5 * [E])).
		KEH_->Mult(hOld, auxRHS);
		KEE_->Mult(eOld, auxFluxdE);
		auxRHS.Add(1.0, auxFluxdE);
		MInv_->Mult(auxRHS, eNew);

		// Update H. dH/dt = M^{-1} * (-S * E + nx * {E} + 0.5 * [H])).
		KHE_->Mult(eOld, auxRHS);
		KHH_->Mult(hOld, auxFluxdH);
		auxRHS.Add(1.0, auxFluxdH);
		MInv_->Mult(auxRHS, hNew);

		break;

	case Centered:
	
		// Update E. dE/dt = M^{-1} * (S * H - {H}).
		KEH_->Mult(hOld, auxRHS);
		MInv_->Mult(auxRHS, eNew);

		// Update H. dH/dt = M^{-1} * (S * E - {E}).
		KHE_->Mult(eOld, auxRHS);
		MInv_->Mult(auxRHS, hNew);

		break;
	}

}

}

