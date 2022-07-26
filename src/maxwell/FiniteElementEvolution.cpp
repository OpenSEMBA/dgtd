#include "FiniteElementEvolution.h"

namespace maxwell {

FiniteElementEvolution::FiniteElementEvolution(FiniteElementSpace* fes, Options options, Model& model, Sources& sources) :
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
			MNND_[f][f2] = buildByMult(buildInverseMassMatrix(f).get(), buildNormalFluxOperator(f2,std::vector<Direction>{}).get());
			for (int dir = Direction::X; dir <= Direction::Z; dir++) {
				Direction d = static_cast<Direction>(dir);
				MS_[f][d] = buildByMult(buildInverseMassMatrix(f).get(), buildDerivativeOperator(d).get());
				MNOD_[f][f2][d] = buildByMult(buildInverseMassMatrix(f).get(), buildNormalFluxOperator(f2, std::vector<Direction>{d}).get());
				if (opts_.disForm == DisForm::Weak) {
					MF_[f][f2][d] = buildByMult(buildInverseMassMatrix(f).get(), buildFluxOperator(opts_.disForm, f2, d).get());
					MP_[f][f2][d] = buildByMult(buildInverseMassMatrix(f).get(), buildPenaltyOperator(opts_.disForm, f2, d).get());
				}
				for (int dir2 = Direction::X; dir2 <= Direction::Z; dir2++) {
					Direction d2 = static_cast<Direction>(dir2);
					MNTD_[f][f2][d][d2] = buildByMult(buildInverseMassMatrix(f).get(), buildNormalFluxOperator(f2, std::vector<Direction>{d, d2}).get());
				}
			}
		}
	}
}

Vector
	FiniteElementEvolution::buildPieceWiseArgVector(const FieldType& f) const
{
	Vector res;
	res.SetSize((int) model_.getAttToMat().size());

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
	FiniteElementEvolution::buildNVector(const Direction& d) const
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

FiniteElementEvolution::FiniteElementOperator
	FiniteElementEvolution::buildInverseMassMatrix(const FieldType& f) const
{
	Vector aux = buildPieceWiseArgVector(f);
	PWConstCoefficient PWCoeff(aux);

	auto MInv = std::make_unique<BilinearForm>(fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));
	MInv->Assemble();
	MInv->Finalize();

	return MInv;
}

FiniteElementEvolution::FiniteElementOperator
	FiniteElementEvolution::buildDerivativeOperator(const Direction& d) const
{
	ConstantCoefficient coeff(0.0);
	auto dir = d;
	auto res = std::make_unique<BilinearForm>(fes_);

	if (d >= fes_->GetMesh()->Dimension()) {
		dir = X;
	}

	switch (opts_.disForm) {
	case DisForm::Weak:
		coeff = ConstantCoefficient(1.0);
		res->AddDomainIntegrator(
			new TransposeIntegrator(
				new DerivativeIntegrator(coeff, dir)
			)
		);
		break;
	case DisForm::Strong: //Hesthaven bilinear form (lambda * u, dvdx). MFEM bilinear form (lambda * dudx, v)
		coeff = ConstantCoefficient(0.5);
		res->AddDomainIntegrator(new DerivativeIntegrator(coeff, dir));
		break;
	}

	res->Assemble();
	res->Finalize();

	if (d >= fes_->GetMesh()->Dimension()) {
		res.get()->operator=(0.0);
	}

	return res;
}

FiniteElementEvolution::FiniteElementOperator
	FiniteElementEvolution::buildFluxOperator(const DisForm& form, const FieldType& f, const Direction& d) const
{
	Vector aux = buildNVector(d);
	VectorConstantCoefficient n(aux);
	auto res = std::make_unique<BilinearForm>(fes_);
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(form, n, c.alpha, c.beta));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(model_.getConstMesh().bdr_attributes.Max());
	for (auto const& kv : model_.getAttToBdr()) {
		Array<int> bdrMarker(model_.getConstMesh().bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int) kv.first - 1] = 1;

		bdrMarkers[(int) kv.first - 1] = bdrMarker;
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.second);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(form, n, c.alpha, c.beta), bdrMarkers[kv.first - 1]);
	}

	res->Assemble();
	res->Finalize();

	if (d >= fes_->GetMesh()->Dimension()) {
		res.get()->operator=(0.0);
	}

	return res;
}

FiniteElementEvolution::FiniteElementOperator
	FiniteElementEvolution::buildNormalFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms) const
{
	std::vector<Direction> dirs = dirTerms;
	auto res = std::make_unique<BilinearForm>(fes_);
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirs, c.beta));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(model_.getConstMesh().bdr_attributes.Max());
	for (auto const& kv : model_.getAttToBdr()) {
		Array<int> bdrMarker(model_.getConstMesh().bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int)kv.first - 1] = 1;

		bdrMarkers[(int)kv.first - 1] = bdrMarker;
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.second);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirs, c.beta), bdrMarkers[kv.first - 1]);
	}

	res->Assemble();
	res->Finalize();

	return res;
}

FiniteElementEvolution::FiniteElementOperator
	FiniteElementEvolution::buildPenaltyOperator(const DisForm& form, const FieldType& f, const Direction& d) const
{

	auto aux = buildNVector(d);
	VectorConstantCoefficient n(aux);
	std::unique_ptr<BilinearForm> res = std::make_unique<BilinearForm>(fes_);
	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(form, n, c.alpha, c.beta));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(model_.getConstMesh().bdr_attributes.Max());
	for (auto const& kv : model_.getAttToBdr()) {
		Array<int> bdrMarker(model_.getConstMesh().bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int) kv.first - 1] = 1;

		bdrMarkers[(int) kv.first - 1] = bdrMarker;
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f, kv.second);
		res->AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(form, n, c.alpha, c.beta), bdrMarkers[kv.first - 1]);
	}

	res->Assemble();
	res->Finalize();

	if (d >= fes_->GetMesh()->Dimension()) {
		res.get()->operator=(0.0);
	}

	return res;
}

FiniteElementEvolution::FiniteElementOperator
	FiniteElementEvolution::buildByMult(
		const BilinearForm* op1, const BilinearForm* op2) const
{
	auto aux = mfem::Mult(op1->SpMat(), op2->SpMat());
	auto res = std::make_unique<BilinearForm>(fes_);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

FiniteElementEvolution::FluxCoefficient
FiniteElementEvolution::interiorFluxCoefficient() const
{
	switch (opts_.disForm) {
	case DisForm::Weak:
		return FluxCoefficient{ 1.0, 0.0 };
	case DisForm::Strong:
		return FluxCoefficient{ 0.0, -0.5 };
	default:
		throw std::exception("No defined BdrCond.");
	}
}
FiniteElementEvolution::FluxCoefficient
	FiniteElementEvolution::interiorPenaltyFluxCoefficient() const
{
	switch (opts_.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		return FluxCoefficient{ 0.0, 0.5 };
	default:
		throw std::exception("No defined FluxType.");
	}
}

FiniteElementEvolution::FluxCoefficient
	FiniteElementEvolution::boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC) const
{
	switch (opts_.disForm) {
	case DisForm::Weak:
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
		default:
			throw std::exception("No defined BdrCond.");
		}
		break;
	case DisForm::Strong:
		switch (bdrC) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0, -2.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0,  0.0 };
			}
		case BdrCond::PMC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0,  0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0, -2.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0, -1.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0, -1.0 };
			}
		default:
			throw std::exception("No defined BdrCond.");
		}
		break;
	default:
		throw std::exception("No defined DisForm.");
	}
}

FiniteElementEvolution::FluxCoefficient
	FiniteElementEvolution::boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC) const
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
				return FluxCoefficient{ 0.0, 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0, 0.0 };
			}
		default:
			throw std::exception("No defined BdrCond.");
		}
	default:
		throw std::exception("No defined FluxType.");
	}

}

void FiniteElementEvolution::Mult(const Vector& in, Vector& out) const
{

	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() + d * fes_->GetNDofs(), fes_->GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_->GetNDofs(), fes_->GetNDofs());
		eNew[d].MakeRef(fes_, &out[d * fes_->GetNDofs()]);
		hNew[d].MakeRef(fes_, &out[(d + 3) * fes_->GetNDofs()]);
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		switch (opts_.disForm) {
		case DisForm::Weak:

			// dtE_x = MS_y * H_z - MF_y * {H_z} - MP_E * [E_z] +
			//        -MS_z * H_y + MF_z * {H_y} + MP_E * [E_y]
			// 
			// Update E.
			MS_[E][z]   ->Mult   (hOld[y], eNew[x]);
			MF_[E][H][z]->AddMult(hOld[y], eNew[x], -1.0);
			MP_[E][E][z]->AddMult(eOld[y], eNew[x], -1.0);
			MS_[E][y]   ->AddMult(hOld[z], eNew[x], -1.0);
			MF_[E][H][y]->AddMult(hOld[z], eNew[x],  1.0);
			MP_[E][E][y]->AddMult(eOld[z], eNew[x],  1.0); 

			// Update H.

			MS_[H][y]   ->Mult   (eOld[z], hNew[x]);
			MF_[H][E][y]->AddMult(eOld[z], hNew[x], -1.0);
			MP_[H][H][y]->AddMult(hOld[z], hNew[x], -1.0);
			MS_[H][z]   ->AddMult(eOld[y], hNew[x], -1.0);
			MF_[H][E][z]->AddMult(eOld[y], hNew[x],  1.0);
			MP_[H][H][z]->AddMult(hOld[y], hNew[x],  1.0);

			break;

		case DisForm::Strong:

			//dtE_x =   MS_y     * Hz   - MS_z * Hy + 
			//		  + MNOD_y   * [Hz] - MNOD_z   * [Hy] + MNND_    * [Ex] +
			//        - MNTD_x_x * [Ex] - MNTD_y_x * [Ey] - MNTD_z_x * [Ez]       
			// Update E.

			MS_[E][y]        ->Mult   (hOld[z], eNew[x]);
			MS_[E][z]        ->AddMult(hOld[y], eNew[x], -1.0);
			MNOD_[E][H][y]   ->AddMult(hOld[z], eNew[x],  1.0);
			MNOD_[E][H][z]   ->AddMult(hOld[y], eNew[x], -1.0);
			MNND_[E][E]		 ->AddMult(eOld[x], eNew[x],  1.0);
			MNTD_[E][E][x][x]->AddMult(eOld[x], eNew[x], -1.0);
			MNTD_[E][E][y][x]->AddMult(eOld[y], eNew[x], -1.0);
			MNTD_[E][E][z][x]->AddMult(eOld[z], eNew[x], -1.0);

			//dtH_x =   MS_z     * Ey   - MS_y     * Ez + 
			//		  + MNOD_z   * [Ey] - MNOD_y   * [Ez] + MNND_    * [Hx] +
			//        - MNTD_x_x * [Hx] - MNTD_y_x * [Hy] - MNTD_z_x * [Hz]    
			//Update H.

			MS_[H][z]		 ->Mult   (eOld[y], hNew[x]);
			MS_[H][y]		 ->AddMult(eOld[z], hNew[x], -1.0);
			MNOD_[H][E][z]	 ->AddMult(eOld[y], hNew[x],  1.0);
			MNOD_[H][E][y]	 ->AddMult(eOld[z], hNew[x], -1.0);
			MNND_[H][H]		 ->AddMult(hOld[x], hNew[x],  1.0);
			MNTD_[H][H][x][x]->AddMult(hOld[x], hNew[x], -1.0);
			MNTD_[H][H][y][x]->AddMult(hOld[y], hNew[x], -1.0);
			MNTD_[H][H][z][x]->AddMult(hOld[z], hNew[x], -1.0);

			break;
		}
	}

}

}

