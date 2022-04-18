#pragma once

#include "mfem.hpp"
#include "BilinearIntegrators.h"

namespace Maxwell1D {

using namespace mfem;

typedef std::size_t Direction;
typedef std::size_t FluxType;
typedef std::size_t Factor;
typedef std::size_t BdrCond;

const Direction X = 0;

const FluxType Centered = 0; //Weak form
const FluxType Upwind = 1; //Strong form

const Factor Alpha = 0;
const Factor Beta = 1;

const BdrCond PEC = 0;
const BdrCond PMC = 1;
const BdrCond SMA = 2;

class FE_Evolution : public TimeDependentOperator {
public:

	static const std::size_t numberOfFieldComponents = 2;
	FluxType fluxType = Upwind;
	BdrCond bdrCond = PEC;

	FE_Evolution(FiniteElementSpace* fes);
	virtual void Mult(const Vector& x, Vector& y) const;
	virtual ~FE_Evolution() = default;

private:

	FiniteElementSpace* fes_;
	
	std::unique_ptr<BilinearForm> MInv_;
	std::unique_ptr<BilinearForm> KEE_; //K - Time Diff. (E) - Applied on (E)
	std::unique_ptr<BilinearForm> KEH_;
	std::unique_ptr<BilinearForm> KHE_;
	std::unique_ptr<BilinearForm> KHH_;

	std::unique_ptr<BilinearForm> buildInverseMassMatrix() const;
	void addDerivativeOperator(
		std::unique_ptr<BilinearForm>& form,
		const Direction& d, ConstantCoefficient& coeff) const;
	void addFluxOperator(
		std::unique_ptr<BilinearForm>& form,
		const Direction& d,const Vector& abgFace,const Vector& abgBdr) const;
	void FE_Evolution::initializeBilinearForms();

};


}