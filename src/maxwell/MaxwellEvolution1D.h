#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "Types.h"
#include "Model.h"
#include "Sources.h"
#include "MaxwellDefs.h"

namespace maxwell {

class MaxwellEvolution1D : public TimeDependentOperator {
public:
	using Vector = mfem::Vector;
	using FiniteElementSpace = mfem::FiniteElementSpace;
	using FiniteElementOperator = mfem::BilinearForm;

	static const int numberOfFieldComponents = 2;

	MaxwellEvolution1D(FiniteElementSpace&, Model&, MaxwellEvolOptions&);
	virtual void Mult(const Vector& x, Vector& y) const;

private:
	FiniteElementSpace& fes_;
	Model& model_;
	MaxwellEvolOptions& opts_;
	
	std::array<FiniteElementOperator, 2> MS_;
	std::array<FiniteElementOperator, 2> MF_;
	std::array<FiniteElementOperator, 2> MP_;
	
};

}