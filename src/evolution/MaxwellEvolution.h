#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "MaxwellEvolutionMethods.h"
#include "components/SubMesher.h"

#include "solver/SourcesManager.h"

#include "components/ProblemDescription.h"
#include "components/DGOperatorFactory.h"

#include <chrono>

namespace maxwell {

void buildGlobalToLocalDoFMapping(const Model& model, const ParFiniteElementSpace& pfes);

class MaxwellEvolution: public mfem::TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	MaxwellEvolution(ProblemDescription& pd, mfem::ParFiniteElementSpace& fes, SourcesManager& srcmngr);
	virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

private:

	std::array<std::unique_ptr<ParBilinearForm>, 2> MInv_;
	
	std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 2> MS_;
	std::array<std::array<std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 3>, 2>, 2> MFNN_;
	std::array<std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 2>, 2> MFN_;
	std::array<std::unique_ptr<ParBilinearForm>, 2> MP_;
	
	std::array<std::unique_ptr<ParBilinearForm>, 2> MPB_;
	std::array<std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 2>, 2> MFNB_;
	std::array<std::array<std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 3>, 2>, 2> MFNNB_;
	
	Vector CND_;
	
	//Total Field and Scattered Field operators for SubMeshing
	std::array<std::unique_ptr<BilinearForm>, 2> MInvTFSF_;
	std::array<std::array<std::unique_ptr<BilinearForm>, 3>, 2> MS_GTFSF_;
	std::array<std::array<std::array<std::array<std::unique_ptr<BilinearForm>, 3>, 3>, 2>, 2> MFNN_GTFSF_;
	std::array<std::array<std::array<std::unique_ptr<BilinearForm>, 3>, 2>, 2> MFN_GTFSF_;
	std::array<std::unique_ptr<BilinearForm>, 2> MP_GTFSF_;

	/* */

	ProblemDescription pd_;
	ParFiniteElementSpace& fes_;
	SourcesManager& srcmngr_;

};

}