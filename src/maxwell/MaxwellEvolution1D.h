#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "Types.h"
#include "Model.h"
#include "Sources.h"
#include "MaxwellDefs.h"

namespace maxwell {

class MaxwellEvolution1D : public TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;

	MaxwellEvolution1D(FiniteElementSpace&, Model&, MaxwellEvolOptions&);
	virtual void Mult(const Vector& x, Vector& y) const;

	const FiniteElementSpace& getFES() { return fes_; }
	std::array<FiniteElementOperator, 2>& getMS() { return MS_; }
	std::array<FiniteElementOperator, 2>& getMF() { return MF_; }
	std::array<FiniteElementOperator, 2>& getMP() { return MP_; }


private:

	std::array<FiniteElementOperator, 2> MS_;
	std::array<FiniteElementOperator, 2> MF_;
	std::array<FiniteElementOperator, 2> MP_;

	FiniteElementSpace& fes_;
	Model& model_;
	MaxwellEvolOptions& opts_;

	
};

}