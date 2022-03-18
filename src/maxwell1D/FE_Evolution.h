#pragma once

#include "mfem.hpp"

namespace Maxwell1D {

using namespace mfem;

typedef std::size_t Direction;
typedef std::size_t FieldType;

const Direction X = 0;

const FieldType Electric = 0;
const FieldType Magnetic = 1;

class FE_Evolution : public TimeDependentOperator {
public:
	static const std::size_t numberOfFieldComponents = 2;

	FE_Evolution(FiniteElementSpace* fes);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:

	FiniteElementSpace* fes_;
	
	std::unique_ptr<BilinearForm> MInv_;
	std::unique_ptr<BilinearForm> KxE_;
	std::unique_ptr<BilinearForm> KxH_;
	std::unique_ptr<BilinearForm> SxE_;
	std::unique_ptr<BilinearForm> SxH_;
	std::unique_ptr<BilinearForm> FxE_;
	std::unique_ptr<BilinearForm> FxH_;

	std::unique_ptr<BilinearForm> buildInverseMassMatrix() const;
	std::unique_ptr<BilinearForm> buildDerivativeAndFluxOperator(
		const Direction& d, const FieldType& ft) const;
	std::unique_ptr<BilinearForm> buildDerivativeOperator(
		const Direction& d, const FieldType& ft) const;
	std::unique_ptr<BilinearForm> buildFluxOperator(
		const Direction& d, const FieldType& ft) const;

};

}