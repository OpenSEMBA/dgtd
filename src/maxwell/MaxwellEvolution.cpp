#include "MaxwellEvolution.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution::FiniteElementOperator
MaxwellEvolution::buildByMult(
	const BilinearForm& op1, 
	const BilinearForm& op2) const
{
	auto aux = mfem::Mult(op1.SpMat(), op2.SpMat());
	auto res = std::make_unique<BilinearForm>(&fes_);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

FieldType altField(const FieldType& f) 
{
	switch (f) {
	case E:
		return H;
	case H:
		return E;
	}
	throw std::runtime_error("Invalid field type for altField.");
}

MaxwellEvolution::MaxwellEvolution(
	FiniteElementSpace& fes, Model& model, const Options& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{model},
	opts_{ options }
{
	for (auto d: {X, Y, Z}) {
		for (auto f : {E, H}) {
			const auto f2{ altField(f) };
			MS_[f][d] = buildByMult(*buildInverseMassMatrix(f), *buildDerivativeOperator(d));
			MF_[f][d] = buildByMult(*buildInverseMassMatrix(f), *buildFluxOperator(f2, d, false));
			MP_[f][d] = buildByMult(*buildInverseMassMatrix(f), *buildFluxOperator(f2, d, true));
		}
	}
}

Vector MaxwellEvolution::buildNVector(const Direction& d) const
{
	const auto dim{ fes_.GetMesh()->Dimension() };
	assert(d < dim);
	
	Vector r(dim);
	r = 0.0;
	r[d] = 1.0;
	return r;
}

MaxwellEvolution::FiniteElementOperator
MaxwellEvolution::buildInverseMassMatrix(const FieldType& f) const
{
	Vector aux{ model_.buildPiecewiseArgVector(f) };
	PWConstCoefficient PWCoeff(aux);

	auto MInv = std::make_unique<BilinearForm>(&fes_);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));

	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}

MaxwellEvolution::FiniteElementOperator
MaxwellEvolution::buildDerivativeOperator(const Direction& d) const
{
	auto res = std::make_unique<BilinearForm>(&fes_);
	
	if (d >= fes_.GetMesh()->Dimension()) {
		res->Assemble();
		res->Finalize();
		return res;
	}

	ConstantCoefficient coeff(1.0);
	res->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(coeff, d)
		)
	);

	res->Assemble();
	res->Finalize();
	return res;
}

MaxwellEvolution::FiniteElementOperator
MaxwellEvolution::buildFluxOperator(const FieldType& f, const Direction& d, bool usePenaltyCoefficients) const
{
	auto res = std::make_unique<BilinearForm>(&fes_);
	if ( d >= fes_.GetMesh()->Dimension()) {
		res->Assemble();
		res->Finalize();
		return res;
	}
	
	Vector aux = buildNVector(d);
	VectorConstantCoefficient n(aux);
	{
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = interiorPenaltyFluxCoefficient();
		}
		else {
			c = interiorFluxCoefficient();
		}
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n, c.alpha, c.beta));
	}

	for (auto& kv : model_.getBoundaryToMarker()) {
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = boundaryPenaltyFluxCoefficient(f, kv.first);
		}
		else {
			c = boundaryFluxCoefficient(f, kv.first);
		}
		res->AddBdrFaceIntegrator(
			new MaxwellDGTraceIntegrator(n, c.alpha, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FluxCoefficient MaxwellEvolution::interiorFluxCoefficient() const
{
	return FluxCoefficient{ 1.0, 0.0 };
}

FluxCoefficient MaxwellEvolution::interiorPenaltyFluxCoefficient() const
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

FluxCoefficient MaxwellEvolution::boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC) const
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
	default:
		throw std::exception("No defined BdrCond.");
	}
}

FluxCoefficient MaxwellEvolution::boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC) const
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
				return FluxCoefficient{ -1.0, 0.0 };
			case FieldType::H:
				return FluxCoefficient{ -1.0, 0.0 };
			}
		default:
			throw std::exception("No defined BdrCond.");
		}
	default:
		throw std::exception("No defined FluxType.");
	}
}

void MaxwellEvolution::Mult(const Vector& in, Vector& out) const
{
	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() + d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		// dtE_x = MS_y * H_z - MF_y * {H_z} - MP_E * [E_z] +
		//        -MS_z * H_y + MF_z * {H_y} + MP_E * [E_y]
		// Update E.
		MS_[E][z]->Mult   (hOld[y], eNew[x]);
		MF_[E][z]->AddMult(hOld[y], eNew[x], -1.0);
		MP_[E][z]->AddMult(eOld[y], eNew[x], -1.0);
		MS_[E][y]->AddMult(hOld[z], eNew[x], -1.0);
		MF_[E][y]->AddMult(hOld[z], eNew[x],  1.0);
		MP_[E][y]->AddMult(eOld[z], eNew[x],  1.0); 

		// Update H.
		MS_[H][y]   ->Mult   (eOld[z], hNew[x]);
		MF_[H][y]->AddMult(eOld[z], hNew[x], -1.0);
		MP_[H][y]->AddMult(hOld[z], hNew[x], -1.0);
		MS_[H][z]   ->AddMult(eOld[y], hNew[x], -1.0);
		MF_[H][z]->AddMult(eOld[y], hNew[x],  1.0);
		MP_[H][z]->AddMult(hOld[y], hNew[x],  1.0);
	}

}

}

