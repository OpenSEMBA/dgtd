#include "FiniteElementEvolution.h"

namespace maxwell {

FiniteElementEvolutionNoCond::FiniteElementEvolutionNoCond(FiniteElementSpace* fes, Options options, Model& model, Sources& sources) :
	TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes->GetNDofs()),
	opts_(options),
	fes_(fes),
	model_(model),
	sources_(sources)
{
	for (int fInt = FieldType::E; fInt <= FieldType::H; fInt++) {
		FieldType f = static_cast<FieldType>(fInt);
		for (int fInt2 = FieldType::E; fInt2 <= FieldType::H; fInt2++) {
			FieldType f2 = static_cast<FieldType>(fInt2);
			for (int dir = Direction::X; dir <= Direction::Z; dir++) {
				Direction d = static_cast<Direction>(dir);
				MS_[f][d] = buildByMult(buildInverseMassMatrix(f).get(), buildDerivativeOperator(d).get());
				MF_[f][f2][d] = buildByMult(buildInverseMassMatrix(f).get(), buildFluxOperator(f2, d).get());
				MP_[f][f2][d] = buildByMult(buildInverseMassMatrix(f).get(), buildPenaltyOperator(f2, d).get());
			}
		}
	}

	setValuesForOperatorsBasedOnProblemDim();
}

void FiniteElementEvolutionNoCond::setValuesForOperatorsBasedOnProblemDim()
{
	if (fes_->GetMesh()->Dimension() == 1){
		for (int i = 0; i < sources_.getSourcesVector().size(); i++) {
			switch (sources_.getSourcesVector().at(i).getDirection()) {
			case X:
				switch (sources_.getSourcesVector().at(i).getFieldType()) {
				case E:
					MS_[E][X].get()->operator=(0.0); MF_[E][H][X].get()->operator=(0.0); MP_[E][E][X].get()->operator=(0.0); MS_[E][Z].get()->operator=(0.0); MF_[E][H][Z].get()->operator=(0.0); MP_[E][E][Z].get()->operator=(0.0);
					MS_[E][Y].get()->operator=(0.0); MF_[E][H][Y].get()->operator=(0.0); MP_[E][E][Y].get()->operator=(0.0); MS_[E][X].get()->operator=(0.0); MF_[E][H][X].get()->operator=(0.0); MP_[E][E][X].get()->operator=(0.0);
					break;
				case H:
					break;
				}
				break;
			case Y:
				switch (sources_.getSourcesVector().at(i).getFieldType()) {
				case E:
					//MS_[E][Y].get()->operator=(0.0); MF_[E][H][Z].get()->operator=(0.0); MP_[E][E][Z].get()->operator=(0.0); MS_[E][Z].get()->operator=(0.0); MF_[E][H][Y].get()->operator=(0.0); MP_[E][E][Y].get()->operator=(0.0);
					//                                 MF_[E][H][X].get()->operator=(0.0); MP_[E][E][X].get()->operator=(0.0);
					/*MS_[E][Z].get()->operator=(0.0); MF_[E][H][X].get()->operator=(0.0);*/ MP_[E][E][X].get()->operator=(0.0); MP_[E][E][Y].get()->operator=(0.0);

					//MS_[H][Z].get()->operator=(0.0); MF_[H][E][Y].get()->operator=(0.0); MP_[H][H][Y].get()->operator=(0.0); MS_[H][Y].get()->operator=(0.0); MF_[H][E][Z].get()->operator=(0.0); MP_[H][H][Z].get()->operator=(0.0);
					//								                                                                                                           
					/*MF_[H][E][X].get()->operator=(0.0); MP_[H][H][X].get()->operator=(0.0);*/ MP_[H][H][X].get()->operator=(0.0); MP_[H][H][Z].get()->operator=(0.0);
					//MS_[H][Y].get()->operator=(0.0); MF_[H][E][X].get()->operator=(0.0); MP_[H][H][X].get()->operator=(0.0);

					break;
				}
				break;
			case Z:
				break;
			}
		}
	}
}

Vector
FiniteElementEvolutionNoCond::buildPieceWiseArgVector(const FieldType& f) const
{
	Vector res;
	res.SetSize(model_.getAttToMat().size());
	
	std::size_t i = 0;
	for (auto const& kv : model_.getAttToMat()) {
		switch (f) {
		case FieldType::E:
			res[i] = kv.second.getPermittivity();
			break;
		case FieldType::H:
			res[i] = kv.second.getPermeability();
			break;
		}
		i++;
	}

	return res;
}

Vector 
FiniteElementEvolutionNoCond::buildNVector(const Direction& d) const
{
	Vector res(fes_->GetMesh()->Dimension());
	switch (fes_->GetMesh()->Dimension()) {
	case 1:
		switch (d) {
		case X:
			res = Vector({ 1.0 });
			break;
		default:
			res = Vector({ 0.0 });
			break;
		}
		break;
	case 2:
		switch (d) {
		case X:
			res = Vector({ 1.0,0.0 });
			break;
		case Y:
			res = Vector({ 0.0,1.0 });
			break;
		default:
			res = Vector({ 0.0,0.0 });
			break;
		}
		break;
	case 3:
		switch (d) {
		case X:
			res = Vector({ 1.0,0.0,0.0 });
			break;
		case Y:
			res = Vector({ 0.0,1.0,0.0 });
			break;
		case Z:
			res = Vector({ 0.0,0.0,1.0 });
			break;
		}
		break;
	}
	return res;
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildInverseMassMatrix(const FieldType& f) const
{
	Vector aux = buildPieceWiseArgVector(f);
	PWConstCoefficient PWCoeff(aux);

	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));
	MInv->Assemble();
	MInv->Finalize();

	return MInv;
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildDerivativeOperator(const Direction& d) const
{
	ConstantCoefficient coeff(1.0);
	auto dir = d;
	auto res = std::make_unique<BilinearForm>(fes_);
	
	if (d >= fes_->GetMesh()->Dimension()) {
		dir = X;
	}

	res->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(coeff, dir)
		)
	);

	res->Assemble();
	res->Finalize();

	if (d >= fes_->GetMesh()->Dimension()) {
		res.get()->operator=(0.0);
	}

	return res;
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildFluxOperator(const FieldType& f, const Direction& d) const
{
	Vector aux = buildNVector(d);
	VectorConstantCoefficient n(aux);
	auto res = std::make_unique<BilinearForm>(fes_);
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	
	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(model_.getConstMesh().bdr_attributes.Max());
	for (auto const& kv : model_.getAttToBdr()) {
		Array<int> bdrMarker(model_.getConstMesh().bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[kv.first-1] = 1;

		bdrMarkers[kv.first-1] = bdrMarker;
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.second);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta), bdrMarkers[kv.first-1]);
	}

	res->Assemble();
	res->Finalize();

	return res;
}

FiniteElementEvolutionNoCond::Operator 
FiniteElementEvolutionNoCond::buildPenaltyOperator(const FieldType& f, const Direction& d) const
{
	
	auto aux = buildNVector(d);
	VectorConstantCoefficient n(aux);
	std::unique_ptr<BilinearForm> res = std::make_unique<BilinearForm>(fes_);
	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}
	
	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(model_.getConstMesh().bdr_attributes.Max());
	for (auto const& kv : model_.getAttToBdr()) {
		Array<int> bdrMarker(model_.getConstMesh().bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[kv.first - 1] = 1;

		bdrMarkers[kv.first - 1] = bdrMarker;
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f, kv.second);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta), bdrMarkers[kv.first - 1]);
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
FiniteElementEvolutionNoCond::boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC) const
{
	switch (bdrC) {
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
FiniteElementEvolutionNoCond::boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC) const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		switch (bdrC) {
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

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		// dtE_x = MS_y * H_z - MF_y * {H_z} - MP_E * [E_z] +
		//        -MS_z * H_y + MF_z * {H_y} + MP_E * [E_y]
		// 
		// Update E.
		MS_[E][y]   ->Mult   (hOld[z], eNew[x]);
		MF_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
		MP_[E][E][y]   ->AddMult(eOld[z], eNew[x], -1.0);
		MS_[E][z]   ->AddMult(hOld[y], eNew[x], -1.0);
		MF_[E][H][z]->AddMult(hOld[y], eNew[x],  1.0);
		MP_[E][E][z]   ->AddMult(eOld[y], eNew[x],  1.0);

		eNew[x].Neg();
		// Update H.

		MS_[H][z]   ->Mult   (eOld[y], hNew[x]);
		MF_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
		MP_[H][H][z]   ->AddMult(hOld[y], hNew[x], -1.0);
		MS_[H][y]   ->AddMult(eOld[z], hNew[x], -1.0);
		MF_[H][E][y]->AddMult(eOld[z], hNew[x],  1.0);
		MP_[H][H][y]   ->AddMult(hOld[z], hNew[x],  1.0);

		hNew[x].Neg();
	}
}

}

