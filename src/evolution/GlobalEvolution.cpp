#include "GlobalEvolution.h"

#include <chrono>

namespace maxwell {

GlobalEvolution::GlobalEvolution(
mfem::ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
mfem::TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
fes_{ fes },
model_{ model },
srcmngr_{ srcmngr },
opts_{ options }
{

	fes_.ExchangeFaceNbrData();

	for (auto d = X; d <= Z; d++){
		eOld_[d].SetSpace(&fes_);
		hOld_[d].SetSpace(&fes_);
	}

	globalOperator_ = std::make_unique<mfem::SparseMatrix>(numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs(), numberOfFieldComponents * numberOfMaxDimensions * (fes_.GetNDofs() + fes_.num_face_nbr_dofs));

#ifdef SHOW_TIMER_INFORMATION
	if (Mpi::WorldRank() == 0){
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}
#endif


	Probes probes;
	if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {

		srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));

		auto globalTFSFfes = srcmngr_.getGlobalTFSFSpace();
		auto tfsfMesh = globalTFSFfes->GetMesh();

		Model tfsfModel = Model(*tfsfMesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
		
		ProblemDescription tfsfpd(tfsfModel, probes, srcmngr_.sources, opts_);
		DGOperatorFactory<FiniteElementSpace> tfsfops(tfsfpd, *globalTFSFfes);

		TFSFOperator_ = tfsfops.buildTFSFGlobalOperator();

		auto src_sm = static_cast<mfem::SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
		mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), sub_to_parent_ids_);
	}

	ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
	DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, fes_);

	globalOperator_ = dgops.buildGlobalOperator();

}

const mfem::Vector buildSingleVectorTFSFFunc(const FieldGridFuncs& func)
{
	mfem::Vector res(6 * func[0][0].Size());
	res.UseDevice(true);
	for (auto f : { E, H }) {
		for (auto d : { X, Y, Z }) {
			for (auto v{ 0 }; v < func[f][d].Size(); v++) {
				res[v + (f * 3 + d) * func[f][d].Size()] = func[f][d][v];
			}
		}
	}
	return res;
}

void assertVectorOnDevice(const Vector &v, const std::string &name)
{
	mfem::StopWatch timerTotal, timerMult, timerTFSF;
	timerTotal.Start();
	timerMult.Start();
	globalOperator_->Mult(in, out);
	timerMult.Stop();

	timerTFSF.Start();
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
				source.reset(nullptr);
				break;
			}
		}
	}
	timerTFSF.Stop();
	timerTotal.Stop();

	std::cout << "Serial times - Total: " << timerTotal.RealTime() * 1000 << "ms, Mult: " << timerMult.RealTime() * 1000 << "ms, TFSF: " << timerTFSF.RealTime() * 1000 << "ms\n";

}


}
