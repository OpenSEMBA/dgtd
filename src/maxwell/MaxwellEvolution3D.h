#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "Types.h"
#include "Model.h"
#include "Sources.h"
#include "SourcesManager.h"
#include "MaxwellDefs.h"

namespace maxwell {

class MaxwellEvolution3D: public mfem::TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	MaxwellEvolution3D(mfem::FiniteElementSpace&, Model&, SourcesManager&, MaxwellEvolOptions&);
	virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

private:

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_;
	std::array<FiniteElementOperator, 2> MP_;

	mfem::FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	MaxwellEvolOptions& opts_;
	

};

}