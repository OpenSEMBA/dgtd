#include "EigenvalueEstimator.h"

using ElementID = int;
using ElementAttribute = int;
using TimeStep = double;
using MaxEigenVal = Eigen::dcomplex;

namespace maxwell {

class OptimizationManager {
public:

	OptimizationManager(mfem::FiniteElementSpace&, Model&);

private:

	void calculateTimeStepForElements();
	void assembleElemIdToElemAtt();

	SubMesh assembleElementBasedSubMesh(ElementID&);
	std::pair<FiniteElementSpace, Model> assembleTimeSteppingReqs(ElementID&);

	mfem::FiniteElementSpace& fes_;
	Model& model_;

	double taggerAttribute_ = 999;

	double optimaldt_;
	std::map<ElementID, ElementAttribute> elemIdToAtt_;
	std::map<ElementID, MaxEigenVal> elemIdToMaxEV_;
};




}