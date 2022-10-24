#include "MaxwellDefs1D.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

Vector buildNVector(const Direction& d, const FiniteElementSpace& fes)
{
	const auto dim{ fes.GetMesh()->Dimension() };
	assert(d < dim);

	Vector r(dim);
	r = 0.0;
	r[d] = 1.0;
	return r;
}

FiniteElementOperator buildFluxOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);
	auto vec = buildNVector(dirTerms.at(0), fes);

	VectorConstantCoefficient vecCC(vec);
	{
		FluxCoefficient c = interiorCenteredFluxCoefficient1D();
		res->AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker())
	{
		FluxCoefficient c = boundaryCenteredFluxCoefficient1D(f, kv.first);
		res->AddBdrFaceIntegrator(
			new MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildPenaltyOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);
	VectorConstantCoefficient one(Vector({ 1.0 }));
	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient1D(opts);
		res->AddInteriorFaceIntegrator(new DGTraceIntegrator(one, c.alpha, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker())
	{
		FluxCoefficient c = boundaryPenaltyFluxCoefficient1D(f, kv.first, opts);
		res->AddBdrFaceIntegrator(
			new DGTraceIntegrator(one, c.alpha, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FluxCoefficient interiorCenteredFluxCoefficient1D()
{
	return FluxCoefficient{ 0.0, 1.0 };
}

FluxCoefficient interiorPenaltyFluxCoefficient1D(const MaxwellEvolOptions& opts)
{
	switch (opts.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		return FluxCoefficient{ 0.0, 1.0 };
	default:
		throw std::exception("No defined FluxType.");
	}
}

FluxCoefficient boundaryCenteredFluxCoefficient1D(const FieldType& f, const BdrCond& bdrC)
{
	switch (bdrC) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0,  0.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0,  2.0 };
		}
	case BdrCond::PMC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0,  2.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0,  0.0 };
		}
	case BdrCond::SMA:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0,  1.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0,  1.0 };
		}
	default:
		throw std::exception("No defined BdrCond.");
	}
}

FluxCoefficient boundaryPenaltyFluxCoefficient1D(const FieldType& f, const BdrCond& bdrC, const MaxwellEvolOptions& opts)
{
	switch (opts.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0, 0.0 };
	case FluxType::Upwind:
		switch (bdrC) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0,  0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0,  2.0 };
			}
		case BdrCond::PMC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0,  2.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0,  0.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0,  1.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0,  1.0 };
			}
		default:
			throw std::exception("No defined BdrCond.");
		}
	}
}
}