#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"

namespace Maxwell1D {

using namespace mfem;

typedef std::size_t Direction;
typedef std::size_t FieldType;
typedef std::size_t FluxType;

const Direction X = 0;

const FieldType Electric = 0;
const FieldType Magnetic = 1;

const FluxType Centered = 0;
const FluxType Upwind = 1;

class FE_Evolution : public TimeDependentOperator {
public:

	static const std::size_t numberOfFieldComponents = 2;
	FluxType fluxType = Upwind;

	FE_Evolution(FiniteElementSpace* fes);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:

	FiniteElementSpace* fes_;
	
	std::unique_ptr<BilinearForm> MInv_;
	std::unique_ptr<BilinearForm> KE_;
	std::unique_ptr<BilinearForm> KH_;
	std::unique_ptr<BilinearForm> SE_;
	std::unique_ptr<BilinearForm> SH_;
	std::unique_ptr<BilinearForm> FE_;
	std::unique_ptr<BilinearForm> FH_;

	std::unique_ptr<BilinearForm> buildInverseMassMatrix() const;
	std::unique_ptr<BilinearForm> buildDerivativeOperator(
		const Direction& d, const FieldType& ft) const;
	std::unique_ptr<BilinearForm> buildFluxOperator(
		const Direction& d, const FieldType& ft) const;

};


}