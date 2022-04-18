#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"

namespace Maxwell1D {

using namespace mfem;

class FE_Evolution : public TimeDependentOperator {
public:
	enum class FluxType {
		Centered,
		Upwind
	};

	enum class BdrCond {
		PEC,
		PMC,
		SMA
	};

	
	struct Options {
		FluxType fluxType = FluxType::Upwind;
		BdrCond bdrCond = BdrCond::PEC;
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
		
	enum class Field {
		Electric,
		Magnetic
	};

	FiniteElementSpace* fes_;
	Options opts_;

	std::unique_ptr<BilinearForm> MInv_, K_;
	FluxOperators FE_, FH_;
	
	void constructBilinearForms();
	std::unique_ptr<BilinearForm> buildInverseMassMatrix() const;
	std::unique_ptr<BilinearForm> buildDerivativeOperator() const;
	FluxOperators buildFluxOperators(const Field&) const;
	
	FluxCoefficient interiorFluxCoefficient() const;
	FluxCoefficient interiorAltFluxCoefficient() const;
	FluxCoefficient boundaryFluxCoefficient(const Field&) const;
	FluxCoefficient boundaryAltFluxCoefficient(const Field&) const;
};


}