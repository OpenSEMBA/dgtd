#include "GlobalEvolution.h"

#include <chrono>

namespace maxwell {

GlobalEvolution::GlobalEvolution(
	ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
	TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options },
	TFSFOperator_{ },
	globalOperator_{ }
{
#ifdef SHOW_TIMER_INFORMATION
	auto startTime{ std::chrono::high_resolution_clock::now() };
#endif

	globalOperator_ = std::make_unique<SparseMatrix>(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(), numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
#endif

	Probes probes;
	if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {

		srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));

		auto globalTFSFfes = srcmngr_.getGlobalTFSFSpace();
		auto tfsfMesh = globalTFSFfes->GetMesh();

		Model tfsfModel = Model(*tfsfMesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
		
		ProblemDescription tfsfpd(tfsfModel, probes, srcmngr_.sources, opts_);
		DGOperatorFactory tfsfops(tfsfpd, *globalTFSFfes);

		TFSFOperator_ = tfsfops.buildTFSFGlobalOperator();

		auto src_sm = static_cast<SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
		SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), sub_to_parent_ids_);
	}

	ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
	DGOperatorFactory dgops(pd, fes_);

	globalOperator_ = dgops.buildGlobalOperator();

}

const Vector buildSingleVectorTFSFFunc(const FieldGridFuncs& func)
{
	Vector res(6 * func[0][0].Size());
	for (auto f : { E, H }) {
		for (auto d : { X, Y, Z }) {
			for (auto v{ 0 }; v < func[f][d].Size(); v++) {
				res[v + (f * 3 + d) * func[f][d].Size()] = func[f][d][v];
			}
		}
	}
	return res;
}

void GlobalEvolution::Mult(const Vector& in, Vector& out) const
{
	const auto& dim{ fes_.GetMesh()->Dimension() };

	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eNew[d].SetSpace(&fes_);
		hNew[d].SetSpace(&fes_);
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
		eNew[d] = 0.0;
		hNew[d] = 0.0;
	}

	globalOperator_->Mult(in, out);

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<TotalField*>(source.get())) {

			auto func{ evalTimeVarFunction(GetTime(),srcmngr_) };
			Vector assembledFunc = buildSingleVectorTFSFFunc(func), tempTFSF(assembledFunc.Size());
			TFSFOperator_->Mult(assembledFunc, tempTFSF);
			for (auto f : { E, H }) {
				for (auto d : { X, Y, Z }) {
					for (auto v{ 0 }; v < sub_to_parent_ids_.Size(); v++) {
						if (tempTFSF[(f * 3 + d) * srcmngr_.getGlobalTFSFSpace()->GetNDofs() + v] != 0.0) {
							out[(f * 3 + d) * fes_.GetNDofs() + sub_to_parent_ids_[v]] -= tempTFSF[(f * 3 + d) * srcmngr_.getGlobalTFSFSpace()->GetNDofs() + v];
						}
					}
				}
			}

		}
	}



}

}