#include "FiniteElementEvolution.h"

namespace maxwell {

//	FiniteElementEvolutionSimple::FiniteElementEvolutionSimple(FiniteElementSpace* fes, Options options) :
//	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
//	opts_(options),
//	fes_(fes),
//	MS_(applyMassOperatorOnOtherOperators(OperatorType::Stiffness)),
//	FEE_(applyMassOperatorOnOtherOperators(OperatorType::Penalty, FieldType::Electric)),
//	FHH_(applyMassOperatorOnOtherOperators(OperatorType::Penalty, FieldType::Magnetic)),
//	FEH_(applyMassOperatorOnOtherOperators(OperatorType::Flux, FieldType::Electric)),
//	FHE_(applyMassOperatorOnOtherOperators(OperatorType::Flux, FieldType::Magnetic))
//{
//}
//
//
//FiniteElementEvolutionSimple::Operator FiniteElementEvolutionSimple::buildInverseMassMatrix() const
//{
//	auto MInv = std::make_unique<BilinearForm>(fes_);
//	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
//	
//	MInv->Assemble();
//	MInv->Finalize();
//	
//	return MInv;
//}
//
//FiniteElementEvolutionSimple::Operator FiniteElementEvolutionSimple::buildDerivativeOperator() const
//{
//	std::size_t d = 0;
//	ConstantCoefficient coeff(1.0);
//
//	auto K = std::make_unique<BilinearForm>(fes_);
//	K->AddDomainIntegrator(
//		new TransposeIntegrator(
//			new DerivativeIntegrator(coeff, d)
//		)
//	);
//
//	K->Assemble();
//	K->Finalize();
//	
//	return K;
//}
//
//FiniteElementEvolutionSimple::Operator FiniteElementEvolutionSimple::buildFluxOperator(const FieldType& f) const
//{
//	auto flux = std::make_unique<BilinearForm>(fes_);
//	VectorConstantCoefficient n(Vector({ 1.0 }));
//	{
//		FluxCoefficient c = interiorFluxCoefficient();
//		flux->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
//	}
//	{
//		FluxCoefficient c = boundaryFluxCoefficient(f);
//		flux->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
//	}
//	flux->Assemble();
//	flux->Finalize();
//
//	return flux;
//}
//
//
//FiniteElementEvolutionSimple::Operator FiniteElementEvolutionSimple::buildPenaltyOperator(const FieldType& f) const
//{
//	std::unique_ptr<BilinearForm> res = std::make_unique<BilinearForm>(fes_);
//
//	VectorConstantCoefficient n(Vector({ 1.0 }));
//	{
//		FluxCoefficient c = interiorPenaltyFluxCoefficient();
//		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
//	}
//	{
//		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f);
//		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
//	}
//
//	res->Assemble();
//	res->Finalize();
//
//	return res;
//}
//
//FiniteElementEvolutionSimple::Operator FiniteElementEvolutionSimple::applyMassOperatorOnOtherOperators(const OperatorType& optype, const FieldType& f) const
//{
//	auto mass = buildInverseMassMatrix();
//	auto second = std::make_unique<BilinearForm>(fes_);
//
//	switch (optype) {
//		case OperatorType::Stiffness:
//			second = buildDerivativeOperator();
//			break;
//		case OperatorType::Flux:
//			second = buildFluxOperator(f);
//			break;
//		case OperatorType::Penalty:
//			second = buildPenaltyOperator(f);
//			break;
//	}
//
//	auto aux = mfem::Mult(mass->SpMat(), second->SpMat());
//
//	auto res = std::make_unique<BilinearForm>(fes_);
//	res->Assemble();
//	res->Finalize();
//	res->SpMat().Swap(*aux);
//
//	return res;
//}
//
//FiniteElementEvolutionSimple::FluxCoefficient FiniteElementEvolutionSimple::interiorFluxCoefficient() const
//{
//	return FluxCoefficient{1.0, 0.0};
//}
//
//FiniteElementEvolutionSimple::FluxCoefficient FiniteElementEvolutionSimple::interiorPenaltyFluxCoefficient() const
//{
//	switch (opts_.fluxType) {
//	case FluxType::Centered:
//		return FluxCoefficient{0.0, 0.0};
//	case FluxType::Upwind:
//		return FluxCoefficient{0.0, 0.5}; 
//	}
//}
//
//FiniteElementEvolutionSimple::FluxCoefficient FiniteElementEvolutionSimple::boundaryFluxCoefficient(const FieldType& f) const
//{
//	switch (opts_.bdrCond) {
//	case BdrCond::PEC:
//		switch (f) {
//		case FieldType::Electric:
//			return FluxCoefficient{ 0.0, 0.0 };
//		case FieldType::Magnetic:
//			return FluxCoefficient{ 2.0, 0.0 };
//		}
//	case BdrCond::PMC:
//		switch (f) {
//		case FieldType::Electric:
//			return FluxCoefficient{ 2.0, 0.0 };
//		case FieldType::Magnetic:
//			return FluxCoefficient{ 0.0, 0.0 };
//		}
//	case BdrCond::SMA:
//		switch (f) {
//		case FieldType::Electric:
//			return FluxCoefficient{ 1.0, 0.0 };
//		case FieldType::Magnetic:
//			return FluxCoefficient{ 1.0, 0.0 };
//		}
//	}
//}
//
//FiniteElementEvolutionSimple::FluxCoefficient FiniteElementEvolutionSimple::boundaryPenaltyFluxCoefficient(const FieldType& f) const
//{
//	switch (opts_.fluxType) {
//	case FluxType::Centered:
//		return FluxCoefficient{ 0.0, 0.0 };
//	case FluxType::Upwind:
//		switch (opts_.bdrCond) {
//		case BdrCond::PEC:
//			switch (f) {
//			case FieldType::Electric:
//				return FluxCoefficient{ 0.0, 0.0 };
//			case FieldType::Magnetic:
//				return FluxCoefficient{ 0.0, 0.0 };
//			}
//		case BdrCond::PMC:
//			switch (f) {
//			case FieldType::Electric:
//				return FluxCoefficient{ 0.0, 0.0 };
//			case FieldType::Magnetic:
//				return FluxCoefficient{ 0.0, 0.0 };
//			}
//		case BdrCond::SMA:
//			switch (f) {
//			case FieldType::Electric:
//				return FluxCoefficient{ 0.0, 0.5 };
//			case FieldType::Magnetic:
//				return FluxCoefficient{ 0.0, 0.5 };
//			}
//		}
//	}
//}
//
//void FiniteElementEvolutionSimple::Mult(const Vector& x, Vector& y) const
//{
//	Vector eOld(x.GetData(), fes_->GetNDofs());
//	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());
//
//	GridFunction eNew(fes_, &y[0]);
//	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);
//
//	// Update E. dE/dt = MS * H - FEH * {H} - FEE * [E].
//
//	MS_->Mult(hOld, eNew);
//	FEH_->AddMult(hOld, eNew, -1.0);
//	FEE_->AddMult(eOld, eNew, -1.0);
//
//	// Update H. dH/dt = MS * E - HE * {E} - FHH * [H].
//
//	MS_->Mult(eOld, hNew);
//	FHE_->AddMult(eOld, hNew, -1.0);
//	FHH_->AddMult(hOld, hNew, -1.0);
//
//}


FiniteElementEvolutionNoCond::FiniteElementEvolutionNoCond(FiniteElementSpace* fes, Options options) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	epsilonVal_(opts_.epsilonVal),
	muVal_(opts_.muVal),
	MES_(
		buildByMult(buildInverseMassMatrix(FieldType::Electric).get(), 
			buildDerivativeOperator(0).get())),
	MHS_(
		buildByMult(buildInverseMassMatrix(FieldType::Magnetic).get(),
			buildDerivativeOperator(0).get())),
	FEE_(buildByMult(buildPenaltyOperator(FieldType::Electric).get(), FieldType::Electric)),
	FHH_(applyMassOperatorOnOtherOperators(buildPenaltyOperator(FieldType::Magnetic).get(), FieldType::Magnetic)),
	FEH_(applyMassOperatorOnOtherOperators(buildFluxOperator(FieldType::Electric).get(), FieldType::Electric)),
	FHE_(applyMassOperatorOnOtherOperators(buildFluxOperator(FieldType::Magnetic).get(), FieldType::Magnetic))

{
	for (int fInt = FieldType::Electric; fInt != FieldType::Magnetic; fInt++) {
		FieldType f = static_cast<FieldType>(fInt);
		for (int d = 0; d < 3; d++) {
			MS_[f][d] = buildByMult(buildInverseMassMatrix(f).get(), buildDerivativeOperator(d).get());
		}
	}
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildInverseMassMatrix(const FieldType& f) const
{
	assert(false); //TODO TODO TODO TODO
	assert(epsilonVal_ != NULL, "epsilonVal Vector is Null");
	Vector aux(epsilonVal_);
	PWConstCoefficient epsilonPWC(aux);

	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(epsilonPWC)));

	MInv->Assemble();
	MInv->Finalize();

	return MInv;
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildDerivativeOperator(Direction d) const
{
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

FiniteElementEvolutionNoCond::Operator FiniteElementEvolutionNoCond::buildFluxOperator(const FieldType& f) const
{

	VectorConstantCoefficient n(Vector({ 1.0 }));
	auto flux = std::make_unique<BilinearForm>(fes_);
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


FiniteElementEvolutionNoCond::Operator FiniteElementEvolutionNoCond::buildPenaltyOperator(const FieldType& f) const
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

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildByMult(
	const BilinearForm* op1, const BilinearForm* op2) const
{
	auto aux = mfem::Mult(op1->SpMat(), op2->SpMat());
	auto res = std::make_unique<BilinearForm>(fes_);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

FiniteElementEvolutionNoCond::FluxCoefficient FiniteElementEvolutionNoCond::interiorFluxCoefficient() const
{
	return FluxCoefficient{ 1.0, 0.0 };
}

FiniteElementEvolutionNoCond::FluxCoefficient FiniteElementEvolutionNoCond::interiorPenaltyFluxCoefficient() const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		return FluxCoefficient{ 0.0, 0.5 };
	}
}

FiniteElementEvolutionNoCond::FluxCoefficient FiniteElementEvolutionNoCond::boundaryFluxCoefficient(const FieldType& f) const
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

FiniteElementEvolutionNoCond::FluxCoefficient FiniteElementEvolutionNoCond::boundaryPenaltyFluxCoefficient(const FieldType& f) const
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

void FiniteElementEvolutionNoCond::Mult(const Vector& x, Vector& y) const
{
	Vector eOld(x.GetData(), fes_->GetNDofs());
	Vector hOld(x.GetData() + fes_->GetNDofs(), fes_->GetNDofs());

	GridFunction eNew(fes_, &y[0]);
	GridFunction hNew(fes_, &y[fes_->GetNDofs()]);

	Vector auxRHS(fes_->GetNDofs());

	// Update E. dE/dt = MS * H - FEH * {H} - FEE * [E]. enew

	MS_->Mult(hOld, eNew);
	FEH_->AddMult(hOld, eNew, -1.0);
	FEE_->AddMult(eOld, eNew, -1.0);

	// Update H. dH/dt = MS * E - HE * {E} - FHH * [H]. hnew

	MS_->Mult(eOld, hNew);
	FHE_->AddMult(eOld, hNew, -1.0);
	FHH_->AddMult(hOld, hNew, -1.0);

}

}

