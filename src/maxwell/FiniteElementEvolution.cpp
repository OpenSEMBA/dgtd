#include "FiniteElementEvolution.h"

namespace maxwell {

FiniteElementEvolutionNoCond::FiniteElementEvolutionNoCond(FiniteElementSpace* fes, Options options, Model& model) :
	TimeDependentOperator(numberOfFieldComponents* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	model_(model)
{
	initializeMaterialParameterVectors();
	getMaterialParameterVectors();

	for (int fInt = FieldType::E; fInt != FieldType::H; fInt++) {
		FieldType f = static_cast<FieldType>(fInt);
		for (int fInt2 = FieldType::E; fInt2 != FieldType::H; fInt2++) {
			FieldType f2 = static_cast<FieldType>(fInt2);
			for (int dir = Direction::X; dir != Direction::Z; dir++) {
				Direction d = static_cast<Direction>(dir);
				MS_[f][d] = buildByMult(buildInverseMassMatrix(f).get(), buildDerivativeOperator(d).get());
				MF_[f][f2][d] = buildByMult(buildInverseMassMatrix(f).get(), buildFluxOperator(f2).get());

			}
			MP_[f][f2] = buildByMult(buildInverseMassMatrix(f).get(), buildPenaltyOperator(f2).get());
		}
	}
}

template<typename KeyType, typename ValueType>
std::pair<KeyType, ValueType> get_max(const std::map<KeyType, ValueType>& x) {
	using pairtype = std::pair<KeyType, ValueType>;
	return *std::max_element(x.begin(), x.end(), [](const pairtype& p1, const pairtype& p2) {
		return p1.second < p2.second;
		});
}

void FiniteElementEvolutionNoCond::initializeMaterialParameterVectors()
{
	std::map<attribute, Material> matMap = model_.getMaterialMap();
	auto maxVals = get_max(matMap);
	eps_.SetSize(maxVals.first);
	mu_.SetSize(maxVals.first);
}

void FiniteElementEvolutionNoCond::getMaterialParameterVectors()
{
	std::map<attribute, Material> matMap = model_.getMaterialMap();
	for(const auto & it : matMap) {
		eps_[it.first-1] = it.second.getPermittivity();
		mu_[it.first-1] = it.second.getPermeability();
	}
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildInverseMassMatrix(const FieldType& f) const
{
	//TODO TODO TODO TODO
	assert(eps_ != NULL, "epsilonVal Vector is Null");
	Vector aux(eps_);
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

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildFluxOperator(const FieldType& f) const
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


FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildPenaltyOperator(const FieldType& f) const
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

FiniteElementEvolutionNoCond::FluxCoefficient 
FiniteElementEvolutionNoCond::interiorFluxCoefficient() const
{
	return FluxCoefficient{ 1.0, 0.0 };
}

FiniteElementEvolutionNoCond::FluxCoefficient 
FiniteElementEvolutionNoCond::interiorPenaltyFluxCoefficient() const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		return FluxCoefficient{ 0.0, 0.5 };
	}
}

FiniteElementEvolutionNoCond::FluxCoefficient 
FiniteElementEvolutionNoCond::boundaryFluxCoefficient(const FieldType& f) const
{
	switch (opts_.bdrCond) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0, 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 2.0, 0.0 };
		}
	case BdrCond::PMC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 2.0, 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0, 0.0 };
		}
	case BdrCond::SMA:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 1.0, 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 1.0, 0.0 };
		}
	}
}

FiniteElementEvolutionNoCond::FluxCoefficient 
FiniteElementEvolutionNoCond::boundaryPenaltyFluxCoefficient(const FieldType& f) const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		switch (opts_.bdrCond) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0, 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		case BdrCond::PMC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0, 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0, 0.5 };
			case FieldType::H:
				return FluxCoefficient{ 0.0, 0.5 };
			}
		}
	}
}

void FiniteElementEvolutionNoCond::Mult(const Vector& in, Vector& out) const
{

	GridFunction eAux;
	
	std::array<Vector,3> eOld, hOld;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() +     d*fes_->GetNDofs(), fes_->GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d+3)*fes_->GetNDofs(), fes_->GetNDofs());
	}

	//if (source) {
	//	eOld[X].projectCoefficient(source.getFunction(t, X, E));
	//}


	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eNew[d].MakeRef(fes_, &out[    d* fes_->GetNDofs()]);
		hNew[d].MakeRef(fes_, &out[(d+3)* fes_->GetNDofs()]);
	}

	Vector auxRHS(fes_->GetNDofs());

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		// Update E.
		MS_[E][y]   ->Mult   (hOld[z], eNew[x]);
		MF_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
		MP_[E][E]   ->AddMult(eOld[x], eNew[x], -1.0);
		MS_[E][z]   ->AddMult(hOld[y], eNew[x], -1.0);
		MF_[E][H][z]->AddMult(hOld[y], eNew[x],  1.0);
		MP_[E][E]   ->AddMult(eOld[x], eNew[x],  1.0);

		// Update H.

		MS_[H][z]   ->Mult   (eOld[y], hNew[x]);
		MF_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
		MP_[H][H]   ->AddMult(hOld[x], hNew[x], -1.0);
		MS_[H][y]   ->AddMult(eOld[z], hNew[x], -1.0);
		MF_[H][E][y]->AddMult(eOld[z], hNew[x],  1.0);
		MP_[H][H]   ->AddMult(hOld[x], hNew[x],  1.0);
	}
}

}

