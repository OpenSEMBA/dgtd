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
    MemoryType mem_type = v.GetMemory().GetMemoryType();
    if (mem_type != MemoryType::DEVICE)
    {
        mfem::out << "Warning: Vector '" << name << "' latest data is NOT on device!\n";
    }
    else
    {
        mfem::out << "Vector '" << name << "' latest data is on DEVICE.\n";
    }
}

void GlobalEvolution::Mult(const mfem::Vector& in, mfem::Vector& out) const
{
    mfem::StopWatch timerTotal, timerExchange, timerAssembleInNew, timerMult, timerTFSF;

	timerTotal.Start();
	
	const auto ndofs = fes_.GetNDofs();
	const auto nbrDofs = fes_.num_face_nbr_dofs;
	const auto blockSize = ndofs + nbrDofs;

    timerExchange.Start();
	#ifdef SEMBA_DGTD_ENABLE_CUDA
	load_in_to_eh_gpu(in, eOld_, hOld_, ndofs);
	#else
	for (auto d = X; d <= Z; d++){
		for (auto v = 0; v < ndofs; v++){
			eOld_[d][v] = in[d * ndofs + v];
			hOld_[d][v] = in[(d+3) * ndofs + v];
		}
	}
	#endif
	for (auto d = X; d <= Z; d++){
		eOld_[d].ExchangeFaceNbrData();
		hOld_[d].ExchangeFaceNbrData();
	}
    timerExchange.Stop();
    
	timerAssembleInNew.Start();
	Vector inNew(6*blockSize);
	#ifdef SEMBA_DGTD_ENABLE_CUDA
	inNew.UseDevice(true);
	load_eh_to_innew_gpu(in, inNew, ndofs, nbrDofs);
	load_nbr_to_innew_gpu(eOld_, hOld_, inNew, ndofs, nbrDofs);
	#else
	for (int d = X; d <= Z; d++) {
		eOld_[d].SetData(in.GetData() + d * ndofs);
		hOld_[d].SetData(in.GetData() + (d + 3) * ndofs);
		inNew.SetVector(eOld_[d],      d  * (ndofs + nbrDofs));
		inNew.SetVector(hOld_[d], (3 + d) * (ndofs + nbrDofs));
		if (nbrDofs != 0){
			inNew.SetVector(eOld_[d].FaceNbrData(),      d  * (ndofs + nbrDofs) + ndofs);
			inNew.SetVector(hOld_[d].FaceNbrData(), (3 + d) * (ndofs + nbrDofs) + ndofs);
		}
	}
	#endif
	timerAssembleInNew.Stop();

	timerMult.Start();
    globalOperator_->Mult(inNew, out);
	timerMult.Stop();

	timerTFSF.Start();
	bool early_tfsf_deletion = false;
    for (auto& source : srcmngr_.sources) {
        if (dynamic_cast<TotalField*>(source.get()) && srcmngr_.getGlobalTFSFSpace() != nullptr) {
			auto sourcecast = dynamic_cast<TotalField*>(source.get());
			if(opts_.tfsfFinalTime != 0.0 
			&& std::abs(GetTime() - opts_.tfsfFinalTime) <= 1e-5){
				#ifdef SEMBA_DGTD_ENABLE_CUDA
				const auto func = eval_time_var_field_gpu(GetTime(), srcmngr_);
            	mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func);
				#else
				auto func = evalTimeVarFunction(GetTime(), srcmngr_);
            	mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
				#endif
				if(assembledFunc.Norml2() >= 1e-1){
					early_tfsf_deletion = true;
					MFEM_WARNING("TFSF SOURCE DELETED WHEN FIELDS ARE STILL BEING LOADED!");
				}
				source.reset(nullptr);
				break;
			}
			#ifdef SEMBA_DGTD_ENABLE_CUDA
            const auto func = eval_time_var_field_gpu(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func);
			#else
			auto func = evalTimeVarFunction(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
			#endif
            mfem::Vector tempTFSF(assembledFunc.Size());
			tempTFSF.UseDevice(true);
            TFSFOperator_->Mult(assembledFunc, tempTFSF);
			#ifdef SEMBA_DGTD_ENABLE_CUDA
			load_tfsf_into_out_vector_gpu(sub_to_parent_ids_, tempTFSF, out, fes_.GetNDofs(), srcmngr_.getGlobalTFSFSpace()->GetNDofs());
			#else
			for (auto f : { E, H }) {
                for (auto d : { X, Y, Z }) {
                    for (int v = 0; v < sub_to_parent_ids_.Size(); v++) {
                        const int outIdx = (f * 3 + d) * ndofs + sub_to_parent_ids_[v];
                        const int tempIdx = (f * 3 + d) * srcmngr_.getGlobalTFSFSpace()->GetNDofs() + v;
                        if (tempTFSF[tempIdx] != 0.0) {
                            out[outIdx] -= tempTFSF[tempIdx];
                        }
                    }
                }
            }
			#endif
		}
	}
	if (early_tfsf_deletion){
		MFEM_WARNING("TFSF SOURCE DELETED WHEN FIELDS ARE STILL BEING LOADED!");
	}
	timerTFSF.Stop();

    timerTotal.Stop();

	std::cout << "Current time: " << GetTime() << std::endl;
    std::cout << "Rank " << Mpi::WorldRank() << " Mult total: " << timerTotal.RealTime() * 1000 << " ms, exchange: " << timerExchange.RealTime() * 1000 << " ms, assembleIn: " << timerAssembleInNew.RealTime() * 1000 << "ms\n";
    std::cout << "Rank " << Mpi::WorldRank() << " Mult mult: " << timerMult.RealTime() * 1000 << " ms, tfsf: " << timerTFSF.RealTime() * 1000 << "ms\n";
}
}