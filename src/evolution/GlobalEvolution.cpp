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

	for (int d = X; d <= Z; d++) {
		eOld_[d].SetSpace(&fes_);
		hOld_[d].SetSpace(&fes_);
		eOld_[d].ExchangeFaceNbrData();
		hOld_[d].ExchangeFaceNbrData();
	}

	inNew_.UseDevice(true);
	inNew_.SetSize(numberOfFieldComponents * numberOfMaxDimensions * (fes_.GetNDofs() + fes_.num_face_nbr_dofs));

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
	for (auto f : { E, H }) {
		for (auto d : { X, Y, Z }) {
			for (auto v{ 0 }; v < func[f][d].Size(); v++) {
				res[v + (f * 3 + d) * func[f][d].Size()] = func[f][d][v];
			}
		}
	}
	return res;
}

void AssertVectorOnDevice(const Vector &v, const std::string &name)
{
    MemoryType mem_type = v.GetMemory().GetMemoryType();
    if (mem_type != MemoryType::DEVICE)
    {
        mfem::out << "Warning: Vector '" << name << "' latest data is NOT on device!\n";
        // You could throw, or forcibly sync device data if you want:
        // const double *dev_ptr = v.DeviceRead();
    }
    else
    {
        mfem::out << "Vector '" << name << "' latest data is on DEVICE.\n";
    }
}

void GlobalEvolution::Mult(const mfem::Vector& in, mfem::Vector& out) const
{
    mfem::StopWatch timerTotal, timerExchange, timerAssembleInNew, timerLoadOutHost, timerMult, timerTFSF;
    timerTotal.Start();
	
	const auto ndofs = fes_.GetNDofs();
	const auto nbr_dofs = fes_.num_face_nbr_dofs;
	const auto block_size = ndofs + nbr_dofs;

	Vector inmod(in);

    timerExchange.Start();
    for (int d = X; d <= Z; d++) {
        eOld_[d].ExchangeFaceNbrData();
        hOld_[d].ExchangeFaceNbrData();
    }
    timerExchange.Stop();
    
	timerAssembleInNew.Start();
	std::array<Vector, 3> eOld, hOld;
	Vector inNew(6*ndofs);

	for (int d = X; d <= Z; d++) {
		eOld[d].SetSize(fes_.GetNDofs());
		hOld[d].SetSize(fes_.GetNDofs());
    }

	inmod.Read();
	for (int d = X; d <= Z; d++) {
		eOld[d].MakeRef(inmod, d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].MakeRef(inmod, (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
		// inNew_.SetVector(eOld_[d].FaceNbrData(),      d  * (fes_.GetNDofs() + fes_.num_face_nbr_dofs) + fes_.GetNDofs());
		// inNew_.SetVector(hOld_[d].FaceNbrData(), (3 + d) * (fes_.GetNDofs() + fes_.num_face_nbr_dofs) + fes_.GetNDofs());
		auto inrw = inNew.HostWrite();
		auto erw = eOld[d].HostRead();
		auto hrw = hOld[d].HostRead();
		assert(inrw != nullptr);
		assert(erw != nullptr);
		assert(hrw != nullptr);
		for (int v = 0; v < ndofs; v++){
			inrw[d * ndofs + v] = erw[v];
			inrw[(d + 3) * ndofs + v] = hrw[v];
		}
	}

	inNew.UseDevice(true);
	out.UseDevice(true);
	inNew.Read();

	AssertVectorOnDevice(inNew, "inNew");
	AssertVectorOnDevice(out, "out");
	timerAssembleInNew.Stop();

	timerMult.Start();
    globalOperator_->Mult(inNew, out);
	timerMult.Stop();

	timerLoadOutHost.Start();
	out.HostRead();
	timerLoadOutHost.Stop();

	timerTFSF.Start();
    for (const auto& source : srcmngr_.sources) {
        if (dynamic_cast<TotalField*>(source.get()) && srcmngr_.getGlobalTFSFSpace() != nullptr) {
            auto func = evalTimeVarFunction(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
            mfem::Vector tempTFSF(assembledFunc.Size());
            TFSFOperator_->Mult(assembledFunc, tempTFSF);

            for (auto f : { E, H }) {
                for (auto d : { X, Y, Z }) {
                    for (int v = 0; v < sub_to_parent_ids_.Size(); v++) {
                        const int outIdx = (f * 3 + d) * fes_.GetNDofs() + sub_to_parent_ids_[v];
                        const int tempIdx = (f * 3 + d) * srcmngr_.getGlobalTFSFSpace()->GetNDofs() + v;
                        if (tempTFSF[tempIdx] != 0.0) {
                            out[outIdx] -= tempTFSF[tempIdx];
                        }
                    }
                }
            }
        }
    }
	timerTFSF.Stop();

    timerTotal.Stop();

	std::cout << "Current time: " << GetTime() << std::endl;
    std::cout << "Rank " << Mpi::WorldRank() << " Mult total: " << timerTotal.RealTime() * 1000 << " ms, exchange: " << timerExchange.RealTime() * 1000 << " ms, assembleIn: " << timerAssembleInNew.RealTime() * 1000 << "ms\n";
    std::cout << "Rank " << Mpi::WorldRank() << " Mult mult: " << timerMult.RealTime() * 1000 << " ms, loadOutHost: " << timerLoadOutHost.RealTime() * 1000 << " ms, tfsf: " << timerTFSF.RealTime() * 1000 << "ms\n";
}
}