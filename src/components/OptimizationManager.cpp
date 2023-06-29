#include "OptimizationManager.h"

namespace maxwell {

OptimizationManager::OptimizationManager(mfem::FiniteElementSpace& fes, Model& model) : 
	fes_(fes),
	model_(model)
{
}

SubMesh OptimizationManager::assembleElementBasedSubMesh(ElementID& id) {
	Mesh m{ model_.getConstMesh() };
	m.SetAttribute(id, taggerAttribute_);
	Array<int> atts(1); atts[0] = taggerAttribute_;
	auto res(SubMesh::CreateFromDomain(m, atts));
	res.SetAttribute(0, elemIdToAtt_[id]);
	return res;
}

std::pair<FiniteElementSpace, Model> OptimizationManager::assembleTimeSteppingReqs(ElementID& id) {
	auto submesh{ assembleElementBasedSubMesh(id) };
	FiniteElementSpace fes(&submesh, fes_.FEColl());
	Model model(submesh, AttributeToMaterial{}, AttributeToBoundary{}, AttributeToInteriorConditions{});
	std::pair<FiniteElementSpace, Model> pair(fes,model);
	return pair;
}

void OptimizationManager::calculateTimeStepForElements() 
{
	for (int e = 0; e < model_.getConstMesh().GetNE(); ++e) {
		
		auto reqs{ assembleTimeSteppingReqs(e) };
		auto ev{ EigenvalueEstimator(reqs.first, reqs.second , EvolutionOptions{})};
		elemIdToMaxEV_.emplace(e, ev.getElementMatrix().eigenvalues().maxCoeff());

	}
}

void OptimizationManager::assembleElemIdToElemAtt()
{
	for (int e = 0; e < model_.getConstMesh().GetNE(); ++e){
		elemIdToAtt_.emplace(e, model_.getConstMesh().GetAttribute(e));
	}
}

}