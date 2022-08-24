#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "Types.h"
#include "Model.h"
#include "Sources.h"

namespace maxwell {

class MaxwellEvolution: public mfem::TimeDependentOperator {
public:
	using Vector = mfem::Vector;
	using FiniteElementSpace = mfem::FiniteElementSpace;
	using BilinearForm = mfem::BilinearForm;
	using FiniteElementOperator = std::unique_ptr<BilinearForm>;

	struct Options {
		FluxType fluxType{ FluxType::Upwind };
	};

	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	MaxwellEvolution(FiniteElementSpace&, Model&, const Options&);
	virtual void Mult(const Vector& x, Vector& y) const;

private:
	FiniteElementSpace& fes_;
	Model& model_;
	Options opts_;
	
	std::array<std::array<FiniteElementOperator, 3>, 2> MS_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MF_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MP_;

	Vector buildNVector(const Direction& d) const;
	
	FiniteElementOperator buildDerivativeOperator(const Direction&) const;
	FiniteElementOperator buildInverseMassMatrix(const FieldType&) const;
	FiniteElementOperator buildFluxOperator(const FieldType&, const Direction&) const;
	FiniteElementOperator buildPenaltyOperator(const FieldType&, const Direction&) const;

	FiniteElementOperator buildByMult(const BilinearForm&, const BilinearForm&) const;

	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorPenaltyFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&, const BdrCond& bdrC) const;
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&, const BdrCond& bdrC) const;
};

}