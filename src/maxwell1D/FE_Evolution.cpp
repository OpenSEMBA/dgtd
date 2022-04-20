#include "FE_Evolution.h"

namespace maxwell1D {

	FE_Evolution::FE_Evolution(FiniteElementSpace* fes, Options options) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	MS_(applyMassOperatorOnOtherOperators(OperatorType::Stiffness)),
	FEE_(applyMassOperatorOnOtherOperators(OperatorType::Penalty, FieldType::Electric)),
	FHH_(applyMassOperatorOnOtherOperators(OperatorType::Penalty, FieldType::Magnetic)),
	FEH_(applyMassOperatorOnOtherOperators(OperatorType::Flux, FieldType::Electric)),
	FHE_(applyMassOperatorOnOtherOperators(OperatorType::Flux, FieldType::Magnetic))
{
}

	

FE_Evolution::Operator FE_Evolution::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	
	MInv->Assemble();
	MInv->Finalize();
	
	return MInv;
}

FE_Evolution::Operator FE_Evolution::buildDerivativeOperator() const
{
	std::size_t d = 0;
	ConstantCoefficient coeff(1.0);

	auto K = std::make_unique<BilinearForm>(fes_);
	K->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(coeff, d)
		)
	);

	K->Assemble();
	K->Finalize();
	
	return K;
}

FE_Evolution::Operator FE_Evolution::buildFluxOperator(const FieldType& f) const
{
	auto flux = std::make_unique<BilinearForm>(fes_);
	VectorConstantCoefficient n(Vector({ 1.0 }));
	{
		FluxCoefficient c = interiorFluxCoefficient();
		flux->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryFluxCoefficient(f);
		flux->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	flux->Assemble();
	flux->Finalize();

	return flux;
}


FE_Evolution::Operator FE_Evolution::buildPenaltyOperator(const FieldType& f) const
{
	std::unique_ptr<BilinearForm> res = std::make_unique<BilinearForm>(fes_);

	VectorConstantCoefficient n(Vector({ 1.0 }));
	{
		FluxCoefficient c = interiorAltFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryAltFluxCoefficient(f);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}

	res->Assemble();
	res->Finalize();

	return res;
}

FE_Evolution::Operator FE_Evolution::applyMassOperatorOnOtherOperators(const OperatorType& optype, const FieldType& f) const
{
	auto mass = buildInverseMassMatrix();
	auto second = std::make_unique<BilinearForm>(fes_);

	switch (optype) {
		case OperatorType::Stiffness:
			second = buildDerivativeOperator();
			break;
		case OperatorType::Flux:
			second = buildFluxOperator(f);
			break;
		case OperatorType::Penalty:
			second = buildPenaltyOperator(f);
			break;
	}

	auto aux = mfem::Mult(mass->SpMat(), second->SpMat());

	auto res = std::make_unique<BilinearForm>(fes_);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

FE_Evolution::FluxCoefficient FE_Evolution::interiorFluxCoefficient() const
{
	return FluxCoefficient{1.0, 0.0};
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

FE_Evolution::FluxCoefficient FE_Evolution::boundaryFluxCoefficient(const FieldType& f) const
{
	switch (opts_.bdrCond) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::Electric:
			return FluxCoefficient{ 0.0, 0.0 };
		case FieldType::Magnetic:
			return FluxCoefficient{ 2.0, 0.0 };
		}
	case BdrCond::SMA:
		switch (f) {
		case FieldType::Electric:
			return FluxCoefficient{ 1.0, 0.0 };
		case FieldType::Magnetic:
			return FluxCoefficient{ 1.0, 0.0 };
		}
	}
}

FE_Evolution::FluxCoefficient FE_Evolution::boundaryAltFluxCoefficient(const FieldType& f) const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		switch (opts_.bdrCond) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::Electric:
				return FluxCoefficient{ 0.0, 0.0 };
			case FieldType::Magnetic:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::Electric:
				return FluxCoefficient{ 0.0, 0.0 };
			case FieldType::Magnetic:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		}
	}
}


//void FE_Evolution::constructBilinearForms()
//{
//	MInv_ = buildInverseMassMatrix();
//	K_ = buildDerivativeOperator();
//	FE_ = buildFluxOperators(FieldType::Electric);
//	FH_ = buildFluxOperators(FieldType::Magnetic);
//}

void FE_Evolution::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(), fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	// Update E. dE/dt = MS * H - FEH * {H} - FEE * [E].

	MS_->Mult(hOld, eNew);
	FEH_->AddMult(hOld, eNew, -1.0);
	FEE_->AddMult(eOld, eNew, -1.0);

	// Update H. dH/dt = MS * E - HE * {E} - FHH * [H].

	MS_->Mult(eOld, hNew);
	FHE_->AddMult(eOld, hNew, -1.0);
	FHH_->AddMult(hOld, hNew, -1.0);

}

}

