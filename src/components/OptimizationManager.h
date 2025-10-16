#include "EigenvalueEstimator.h"

using ElementID = int;
using ElementAttribute = int;
using TimeStep = double;
using MaxEigenVal = Eigen::dcomplex;

namespace maxwell {

class OptimizationManager {
public:

	OptimizationManager(mfem::ParFiniteElementSpace&, Model&);

private:

	void calculateTimeStepForElements();
	void assembleElemIdToElemAtt();

	ParSubMesh assembleElementBasedSubMesh(ElementID&);
	std::pair<ParFiniteElementSpace, Model> assembleTimeSteppingReqs(ElementID&);

	mfem::ParFiniteElementSpace& fes_;
	Model& model_;

	int taggerAttribute_ = 999;

	double optimaldt_ = 0.0;
	std::map<ElementID, ElementAttribute> elemIdToAtt_;
	std::map<ElementID, MaxEigenVal> elemIdToMaxEV_;
};




}