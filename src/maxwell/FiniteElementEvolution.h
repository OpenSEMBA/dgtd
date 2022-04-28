#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"
#include "Types.h"

#include <array>

namespace maxwell {

using namespace mfem;

//class FiniteElementEvolutionSimple : public TimeDependentOperator {
//public:
//
//	struct Options {
//		FluxType fluxType = FluxType::Upwind;
//		BdrCond bdrCond = BdrCond::PEC;
//	};
//
//	enum class OperatorType {
//		Stiffness,
//		Flux,
//		Penalty
//	};
//
//	static const std::size_t numberOfFieldComponents = 2;
//	
//	FiniteElementEvolutionSimple(FiniteElementSpace* fes, Options options);
//	virtual void Mult(const Vector& x, Vector& y) const;
//	virtual ~FiniteElementEvolutionSimple() = default;
//
//
//private:
//	struct FluxCoefficient {
//		double alpha;
//		double beta;
//	};
//
//	typedef std::pair<std::unique_ptr<BilinearForm>, std::unique_ptr<BilinearForm>> FluxOperators;
//	typedef std::unique_ptr<BilinearForm> Operator;
//
//	FiniteElementSpace* fes_;
//	Options opts_;
//
//	Operator MS_, FEE_, FEH_, FHE_, FHH_;
//
//	Vector eps_, mu_;
//	const Vector epsilonVal_;
//	const Vector muVal_;
//
//	Operator buildInverseMassMatrix() const;
//	Operator buildDerivativeOperator() const;
//	Operator buildFluxOperator(const FieldType&) const;
//	Operator buildPenaltyOperator(const FieldType&) const;
//
//	Operator applyMassOperatorOnOtherOperators(const OperatorType&, const FieldType& f = FieldType::Electric) const;
//	
//	FluxCoefficient interiorFluxCoefficient() const;
//	FluxCoefficient interiorPenaltyFluxCoefficient() const;
//	FluxCoefficient boundaryFluxCoefficient(const FieldType&) const;
//	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&) const;
//};

class FiniteElementEvolutionNoCond : public TimeDependentOperator {
public:

	struct Options {
		FluxType fluxType = FluxType::Upwind;
		BdrCond bdrCond = BdrCond::PEC;
		Vector epsilonVal = Vector({ 1.0 });
		Vector muVal = Vector({ 1.0 });
	};

	static const std::size_t numberOfFieldComponents = 2;

	FiniteElementEvolutionNoCond(FiniteElementSpace* fes, Options options);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FiniteElementEvolutionNoCond() = default;

private:
	struct FluxCoefficient {
		double alpha;
		double beta;
	};

	typedef std::unique_ptr<BilinearForm> Operator;

	FiniteElementSpace* fes_;
	Options opts_;
	Vector eps_, mu_;
	Vector epsilonVal_;
	Vector muVal_;

	std::array<std::array<Operator, 2>,3> MS_, FE_, FH_;
	
	Operator buildDerivativeOperator(Direction) const;
	Operator buildInverseMassMatrix(const FieldType&) const;
	Operator buildFluxOperator(const FieldType&) const;
	Operator buildPenaltyOperator(const FieldType&) const;

	Operator buildByMult(const BilinearForm*, const BilinearForm*) const;

	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorPenaltyFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&) const;
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&) const;
};

}