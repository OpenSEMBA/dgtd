#include "FiniteElementEvolutionMaterials.h"

namespace maxwell {

	FiniteElementEvolutionMaterials::FiniteElementEvolutionMaterials(FiniteElementSpace* fes, Options options) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	MS_(applyMassOperatorOnOtherOperators(OperatorType::Stiffness)),
	FEE_(applyMassOperatorOnOtherOperators(OperatorType::Penalty, FieldType::Electric)),
	FHH_(applyMassOperatorOnOtherOperators(OperatorType::Penalty, FieldType::Magnetic)),
	FEH_(applyMassOperatorOnOtherOperators(OperatorType::Flux, FieldType::Electric)),
	FHE_(applyMassOperatorOnOtherOperators(OperatorType::Flux, FieldType::Magnetic)),
	epsilonVal_(opts_.epsilonVal),
	muVal_(opts_.muVal)
{
}

Vector FiniteElementEvolutionMaterials::buildEpsilonVector() {
	Vector res(fes_->GetNDofs());
	for (int i = 0; i < fes_->GetMesh()->GetNE(); i++) {
		Array<int> aux;
		fes_->GetElementDofs(i, aux);
		for (int j = 0; j < aux.Size(); j++) {
			res[aux[j]] = epsilonVal_[fes_->GetMesh()->GetAttribute(i) - 1];
		}
	}
	return res;
}

Vector FiniteElementEvolutionMaterials::buildMuVector() {
	Vector res(fes_->GetNDofs());
	for (int i = 0; i < fes_->GetMesh()->GetNE(); i++) {
		Array<int> aux;
		fes_->GetElementDofs(i, aux);
		for (int j = 0; j < aux.Size(); j++) {
			res[aux[j]] = muVal_[fes_->GetMesh()->GetAttribute(i) - 1];
		}
	}
	return res;
}
FiniteElementEvolutionMaterials::Operator FiniteElementEvolutionMaterials::buildInverseMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	
	MInv->Assemble();
	MInv->Finalize();
	
	return MInv;
}

FiniteElementEvolutionMaterials::Operator FiniteElementEvolutionMaterials::buildDerivativeOperator() const
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

FiniteElementEvolutionMaterials::Operator FiniteElementEvolutionMaterials::buildFluxOperator(const FieldType& f) const
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


FiniteElementEvolutionMaterials::Operator FiniteElementEvolutionMaterials::buildPenaltyOperator(const FieldType& f) const
{
	std::unique_ptr<BilinearForm> res = std::make_unique<BilinearForm>(fes_);

	VectorConstantCoefficient n(Vector({ 1.0 }));
	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	{
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}

	res->Assemble();
	res->Finalize();

	return res;
}

FiniteElementEvolutionMaterials::Operator FiniteElementEvolutionMaterials::applyMassOperatorOnOtherOperators(const OperatorType& optype, const FieldType& f) const
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

FiniteElementEvolutionMaterials::FluxCoefficient FiniteElementEvolutionMaterials::interiorFluxCoefficient() const
{
	return FluxCoefficient{1.0, 0.0};
}

FiniteElementEvolutionMaterials::FluxCoefficient FiniteElementEvolutionMaterials::interiorPenaltyFluxCoefficient() const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{0.0, 0.0};
	case FluxType::Upwind:
		return FluxCoefficient{0.0, 0.5}; 
	}
}

FiniteElementEvolutionMaterials::FluxCoefficient FiniteElementEvolutionMaterials::boundaryFluxCoefficient(const FieldType& f) const
{
	switch (opts_.bdrCond) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::Electric:
			return FluxCoefficient{ 0.0, 0.0 };
		case FieldType::Magnetic:
			return FluxCoefficient{ 2.0, 0.0 };
		}
	case BdrCond::PMC:
		switch (f) {
		case FieldType::Electric:
			return FluxCoefficient{ 2.0, 0.0 };
		case FieldType::Magnetic:
			return FluxCoefficient{ 0.0, 0.0 };
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

FiniteElementEvolutionMaterials::FluxCoefficient FiniteElementEvolutionMaterials::boundaryPenaltyFluxCoefficient(const FieldType& f) const
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
		case BdrCond::PMC:
			switch (f) {
			case FieldType::Electric:
				return FluxCoefficient{ 0.0, 0.0 };
			case FieldType::Magnetic:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::Electric:
				return FluxCoefficient{ 0.0, 0.5 };
			case FieldType::Magnetic:
				return FluxCoefficient{ 0.0, 0.5 };
			}
		}
	}
}

void FiniteElementEvolutionMaterials::Mult(const Vector& x, Vector& y) const
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

