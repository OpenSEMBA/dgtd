#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"
#include "Types.h"
#include "Model.h"
#include "Sources.h"

#include <array>

namespace maxwell {

class FiniteElementEvolution : public TimeDependentOperator {
public:

	struct Options {
		FluxType fluxType = FluxType::Upwind;
		DisForm disForm = DisForm::Weak;
	};

	static const std::size_t numberOfFieldComponents = 2;
	static const std::size_t numberOfMaxDimensions = 3;

	FiniteElementEvolution(FiniteElementSpace* fes, Options options, Model& model, Sources& sources);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FiniteElementEvolution() = default;

private:
	struct FluxCoefficient {
		double alpha;
		double beta;
	};

	typedef std::unique_ptr<BilinearForm> Operator;

	FiniteElementSpace* fes_;
	Options opts_;
	Model model_;
	Sources sources_;

	std::array<std::array<Operator, 3>, 2> MS_;
	std::array<std::array<std::array<Operator, 3>, 2>, 2> MF_;
	std::array<std::array<std::array<Operator, 3>, 2>, 2> MP_;

	Vector buildNVector(const Direction& d) const;
	Vector buildPieceWiseArgVector(const FieldType& f) const;
	
	Operator buildDerivativeOperator(const Direction&) const;
	Operator buildInverseMassMatrix(const FieldType&) const;
	Operator buildWeakFluxOperator(const FieldType&, const Direction& d) const;
	Operator buildStrongFluxOperator(const FieldType&, const Direction& d) const;
	Operator buildWeakPenaltyOperator(const FieldType& f, const Direction& d) const;
	Operator buildStrongPenaltyOperator(const FieldType& f, const Direction& d) const;

	Operator buildByMult(const BilinearForm*, const BilinearForm*) const;

	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorPenaltyFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&, const BdrCond& bdrC) const;
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&, const BdrCond& bdrC) const;
};

}