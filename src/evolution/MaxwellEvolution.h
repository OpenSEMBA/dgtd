#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "MaxwellEvolutionMethods.h"
#include "components/SubMesher.h"

#include "solver/SourcesManager.h"

#include "components/ProblemDefinition.h"
#include "components/DGOperatorFactory.h"

#include <chrono>

namespace maxwell {

class MaxwellEvolution: public mfem::TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	MaxwellEvolution(ProblemDescription& pd, mfem::ParFiniteElementSpace& fes, SourcesManager& srcmngr);
	virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

private:

	std::array<FiniteElementOperator, 2> MInv_;
	std::array<FiniteElementOperator, 2> MInvTFSF_;

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_;
	std::array<FiniteElementOperator, 2> MP_;

	std::array<FiniteElementOperator, 2> MPB_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFNB_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNNB_;
	
	Vector CND_;

	//Total Field and Scattered Field operators for SubMeshing

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_GTFSF_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_GTFSF_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_GTFSF_;
	std::array<FiniteElementOperator, 2> MP_GTFSF_;

	/* */

	ProblemDescription pd_;
	FiniteElementSpace& fes_;
	SourcesManager& srcmngr_;

};

}