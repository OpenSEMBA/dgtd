#include "MaxwellEvolutionMethods.h"

namespace maxwell {

using namespace mfem;

FiniteElementOperator buildByMult(
	const BilinearForm& op1,
	const BilinearForm& op2,
	FiniteElementSpace& fes)
{
	auto aux = mfem::Mult(op1.SpMat(), op2.SpMat());
	auto res = std::make_unique<BilinearForm>(&fes);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

Vector buildNVector(const Direction& d, const FiniteElementSpace& fes)
{
	const auto dim{ fes.GetMesh()->Dimension() };
	Vector r(dim);
	r = 0.0;
	if (d < dim) {
		r[d] = 1.0;
		return r;
	}
	else {
		return r;
	}
}

FiniteElementOperator buildInverseMassMatrix(const FieldType& f, const Model& model, FiniteElementSpace& fes)
{
	Vector aux{ model.buildPiecewiseArgVector(f) };
	PWConstCoefficient PWCoeff(aux);

	auto MInv = std::make_unique<BilinearForm>(&fes);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));
		
	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}



FiniteElementOperator buildDerivativeOperator(const Direction& d, FiniteElementSpace& fes)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	if (d >= fes.GetMesh()->Dimension()) {
		res->Assemble();
		res->Finalize();
		return res;
	}

	ConstantCoefficient coeff(1.0);
	res->AddDomainIntegrator(
		new DerivativeIntegrator(coeff, d)
	);

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.first);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildPenaltyOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient(opts);
		res->AddInteriorFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	{
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = interiorPenaltyFluxCoefficient(opts);
		}
		else {
			c = interiorFluxCoefficient();
		}
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		}
		else {
			c = boundaryFluxCoefficient(f, kv.first);
		}
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxJumpOperator(const FieldType& f, const std::vector<Direction>& dirTerms, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	{
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = interiorPenaltyFluxCoefficient(opts);
		}
		else {
			c = interiorFluxCoefficient();
		}
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		}
		else {
			c = boundaryFluxCoefficient(f, kv.first);
		}
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFunctionOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker())
	{
		res->AddInteriorBoundaryFaceIntegrator(
			new mfemExtension::MaxwellDGTotalFieldIntegrator(dirTerms, 1.0), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<BilinearForm>(&fes);
	auto vec = buildNVector(dirTerms.at(0), fes);

	VectorConstantCoefficient vecCC(vec);
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker())
	{
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.first);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second);
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
		FluxCoefficient c = interiorPenaltyFluxCoefficient(opts);
		res->AddInteriorFaceIntegrator(new DGTraceIntegrator(one, 0.0, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker())
	{
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		res->AddBdrFaceIntegrator(
			new DGTraceIntegrator(one, 0.0, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFunctionOperator1D(const FieldType& f, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker())
	{
		res->AddInteriorBoundaryFaceIntegrator(
			new mfemExtension::MaxwellDGTotalFieldIntegrator({ X }, 1.0), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementVector buildBoundaryFunctionVector1D(Model& model, FiniteElementSpace& fes) 
{
	auto res = std::make_unique<mfemExtension::LinearForm>(&fes);
	VectorConstantCoefficient one(Vector({ 1.0 }));

	for (auto& kv : model.getBoundaryToMarker())
	{
		res->AddBdrFaceIntegrator(
			new mfemExtension::BoundaryDGJumpIntegrator(one, 1.0), kv.second
		);
	}

	res->Assemble();
	return res;

}


FluxCoefficient interiorFluxCoefficient()
{
	return FluxCoefficient{ 1.0 };
}

FluxCoefficient interiorPenaltyFluxCoefficient(const MaxwellEvolOptions& opts)
{
	switch (opts.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0 };
	case FluxType::Upwind:
		return FluxCoefficient{ 1.0 };
	default:
		throw std::exception("No defined FluxType.");
	}
}

FluxCoefficient boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC)
{
	switch (bdrC) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 2.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0 };
		}
	case BdrCond::PMC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 2.0 };
		}
	case BdrCond::SMA:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 1.0 };
		case FieldType::H:
			return FluxCoefficient{ 1.0 };
		}
	case BdrCond::TotalField:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0 };
		}	
	case BdrCond::NONE:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0 };
		}
	default:
		throw std::exception("No defined BdrCond.");
	}
}

FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC, const MaxwellEvolOptions& opts)
{
	switch (opts.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0 };
	case FluxType::Upwind:
		switch (bdrC) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 2.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0 };
			}
		case BdrCond::PMC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 2.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 1.0 };
			case FieldType::H:
				return FluxCoefficient{ 1.0 };
			}
		case BdrCond::TotalField:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0 };
			}
		case BdrCond::NONE:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0 };
			}
		default:
			throw std::exception("No defined BdrCond.");
		}
	default:
		throw std::exception("No defined FluxType.");
	}
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

}