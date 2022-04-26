#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"
#include "Types.h"

namespace maxwell {

using namespace mfem;

class FiniteElementEvolutionMaterials : public TimeDependentOperator {
public:

	struct Options {
		FluxType fluxType = FluxType::Upwind;
		BdrCond bdrCond = BdrCond::PEC;
		Vector epsilonVal;
		Vector muVal;
	};

	enum class OperatorType {
		Stiffness,
		Flux,
		Penalty
	};

	static const std::size_t numberOfFieldComponents = 2;

	FiniteElementEvolutionMaterials(FiniteElementSpace* fes, Options options);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FiniteElementEvolutionMaterials() = default;


private:
	struct FluxCoefficient {
		double alpha;
		double beta;
	};

	typedef std::unique_ptr<BilinearForm> Operator;

	FiniteElementSpace* fes_;
	Options opts_;

	Operator MS_, FEE_, FEH_, FHE_, FHH_;

	Vector eps_, mu_;
	const Vector epsilonVal_;
	const Vector muVal_;

	Vector buildEpsilonVector();
	Vector buildMuVector();

	Operator buildInverseMassMatrix() const;
	Operator buildDerivativeOperator() const;
	Operator buildFluxOperator(const FieldType&) const;
	Operator buildPenaltyOperator(const FieldType&) const;

	Operator applyMassOperatorOnOtherOperators(const OperatorType&, const FieldType& f = FieldType::Electric) const;

	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorPenaltyFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&) const;
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&) const;
};


}