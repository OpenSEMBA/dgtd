#include "FE_Evolution.h"

namespace Maxwell1D {

FE_Evolution::FE_Evolution(FiniteElementSpace* fes, Options options) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	MInv_(buildInverseMassMatrix())
{
	constructBilinearForms();
}

void FE_Evolution::assembleAndFinalizeAllOperators()
{
	MInv_->Assemble();
	MInv_->Finalize();
	
	K_->Assemble();
	K_->Finalize();

	FE_.first->Assemble();
	FE_.first->Finalize();
	FE_.second->Assemble();
	FE_.second->Finalize();

	FH_.first->Assemble();
	FH_.first->Finalize();
	FH_.second->Assemble();
	FH_.second->Finalize();
}

std::unique_ptr<BilinearForm> FE_Evolution::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	return MInv;
}

std::unique_ptr<BilinearForm> FE_Evolution::buildDerivativeOperator() const
{
	std::size_t d = 0;
	ConstantCoefficient coeff(1.0);
	auto K = std::make_unique<BilinearForm>(fes_);
	K->AddDomainIntegrator(new TransposeIntegrator(new DerivativeIntegrator(coeff, d)));
	return K;
}

FE_Evolution::FluxOperators FE_Evolution::buildFluxOperators(const Field& f) const
{
	FluxOperators res = std::make_pair(
		std::make_unique<BilinearForm>(fes_),
		std::make_unique<BilinearForm>(fes_));

	VectorConstantCoefficient n = VectorConstantCoefficient(Vector({1.0}));
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res.first->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryFluxCoefficient(f);
		res.first->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = interiorAltFluxCoefficient();
		res.second->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryAltFluxCoefficient(f);
		res.second->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	return res;
}

FE_Evolution::FluxCoefficient FE_Evolution::interiorFluxCoefficient() const
{
	return FluxCoefficient{-1.0, 0.0};
}

FE_Evolution::FluxCoefficient FE_Evolution::interiorAltFluxCoefficient() const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{0.0, 0.0};
	case FluxType::Upwind:
		return FluxCoefficient{0.0, 0.5};
	}
}

FE_Evolution::FluxCoefficient FE_Evolution::boundaryFluxCoefficient(const Field& f) const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0,0.0 };
	case FluxType::Upwind:
		switch (opts_.bdrCond) {
		case BdrCond::PEC:
			switch (f) {
			case Field::Electric:
				return FluxCoefficient{ 0.0, 0.0 };
			case Field::Magnetic:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		}
	}
}

FE_Evolution::FluxCoefficient FE_Evolution::boundaryAltFluxCoefficient(const Field& f) const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0,0.0 };
	case FluxType::Upwind:
		switch (opts_.bdrCond) {
		case BdrCond::PEC:
			switch (f) {
			case Field::Electric:
				return FluxCoefficient{ 0.0, 0.0 };
			case Field::Magnetic:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		}
	}
}


void FE_Evolution::constructBilinearForms()
{
	MInv_ = buildInverseMassMatrix();
	K_ = buildDerivativeOperator();
	FE_ = buildFluxOperators(Field::Electric);
	FH_ = buildFluxOperators(Field::Magnetic);
	assembleAndFinalizeAllOperators();
}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(), fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	Vector auxRHS(MInv_->Height());

	// Update E. dE/dt = M^{-1} * (-K * H + FE * {H} + altFE * [E])).
	auxRHS = 0.0;
	K_->AddMult(hOld, auxRHS, -1.0);
	FE_.first->AddMult(hOld, auxRHS);
	FE_.second->AddMult(eOld, auxRHS);
	MInv_->Mult(auxRHS, eNew);

	// Update H. dH/dt = M^{-1} * (-K * E + FH * {E} + altFH * [H])).
	auxRHS = 0.0;
	K_->AddMult(eOld, auxRHS, -1.0);
	FH_.first->AddMult(eOld, auxRHS);
	FH_.second->AddMult(hOld, auxRHS);
	MInv_->Mult(auxRHS, hNew);

}

}

