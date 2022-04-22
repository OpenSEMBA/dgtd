#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"
#include "Types.h"

namespace maxwell {

using namespace mfem;

class FE_Evolution : public TimeDependentOperator {
public:

	struct Options {
		FluxType fluxType = FluxType::Upwind;
		BdrCond bdrCond = BdrCond::PEC;
	};

	enum class OperatorType {
		Stiffness,
		Flux,
		Penalty
	};
	static const std::size_t numberOfFieldComponents = 2;
	
	FE_Evolution(FiniteElementSpace* fes, Options options);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:
	struct FluxCoefficient {
		double alpha;
		double beta;
	};


	typedef std::pair<std::unique_ptr<BilinearForm>, std::unique_ptr<BilinearForm>> FluxOperators;
	typedef std::unique_ptr<BilinearForm> Operator;

	FiniteElementSpace* fes_;
	Options opts_;

	Operator MS_, FEE_, FEH_, FHE_, FHH_;
	
	void constructBilinearForms();

	Operator buildInverseMassMatrix() const;
	Operator buildDerivativeOperator() const;
	Operator buildFluxOperator(const FieldType&) const;
	Operator buildPenaltyOperator(const FieldType&) const;

	Operator applyMassOperatorOnOtherOperators(const OperatorType&, const FieldType& f = FieldType::Electric) const;

	//Operator buildMassAndStiffOperator() const;
	//Operator buildMassAndFluxOperator(const FieldType&) const;
	//Operator buildMassAndPenaltyOperator(const FieldType&) const;
	
	FluxOperators buildFluxOperators(const FieldType&) const;
	
	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorPenaltyFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const FieldType&) const;
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&) const;
};


}