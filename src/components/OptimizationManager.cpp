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
	Model model(submesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
	std::pair<FiniteElementSpace, Model> pair(fes,model);
	return pair;
}

Eigen::dcomplex calcHighestModulus(Eigen::VectorXcd evs)
{
	Eigen::dcomplex res(0, 0);
	for (int i = 0; i < evs.size(); ++i) {
		if (sqrt(pow(res.real(), 2.0) + pow(res.imag(), 2.0)) < sqrt(pow(evs[i].real(), 2.0) + pow(evs[i].imag(), 2.0))) {
			res = evs[i];
		}
	}
	return res;
}

void OptimizationManager::calculateTimeStepForElements() 
{
	for (int e = 0; e < model_.getConstMesh().GetNE(); ++e) {
		
		auto reqs{ assembleTimeSteppingReqs(e) };
		auto eo{ EvolutionOptions{}};
		auto ev{ EigenvalueEstimator(reqs.first, reqs.second, eo)};
		auto evs{ ev.getElementMatrix().eigenvalues() };
		elemIdToMaxEV_.emplace(e, calcHighestModulus(evs));

	}
}

void OptimizationManager::assembleElemIdToElemAtt()
{
	for (int e = 0; e < model_.getConstMesh().GetNE(); ++e){
		elemIdToAtt_.emplace(e, model_.getConstMesh().GetAttribute(e));
	}
}

}
