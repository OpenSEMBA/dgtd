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
    mfem::StopWatch timerTotal, timerExchange, timerAssembleInNew, timerApplyA, timerTFSF;

    timerTotal.Start();

    const auto ndofs     = fes_.GetNDofs();
    const auto nbrDofs   = fes_.num_face_nbr_dofs;
    const auto blockSize = ndofs + nbrDofs;

    // ---------------------------
    // 1) Load locals and exchange neighbors
    // ---------------------------
    timerExchange.Start();
#ifdef SEMBA_DGTD_ENABLE_CUDA
    load_in_to_eh_gpu(in, eOld_, hOld_, ndofs);
#else
    for (int d = X; d <= Z; ++d)
    {
        for (int v = 0; v < ndofs; ++v)
        {
            eOld_[d][v] = in[d * ndofs + v];
            hOld_[d][v] = in[(3 + d) * ndofs + v];
        }
    }
#endif
    for (int d = X; d <= Z; ++d)
    {
        eOld_[d].ExchangeFaceNbrData();
        hOld_[d].ExchangeFaceNbrData();
    }
    timerExchange.Stop();

    // ---------------------------
    // 2) Assemble packed input [locals|neighbors]
    // ---------------------------
    timerAssembleInNew.Start();
    mfem::Vector inNew(6 * blockSize);
    inNew.UseDevice(true);

#ifdef SEMBA_DGTD_ENABLE_CUDA
    load_eh_to_innew_gpu(in, inNew, ndofs, nbrDofs);
    load_nbr_to_innew_gpu(eOld_, hOld_, inNew, ndofs, nbrDofs);
#else
    for (int d = X; d <= Z; ++d)
    {
        MFEM_ASSERT(eOld_[d].Size() == ndofs, "eOld size mismatch");
        MFEM_ASSERT(hOld_[d].Size() == ndofs, "hOld size mismatch");

        // locals
        inNew.SetVector(eOld_[d],       d      * blockSize);
        inNew.SetVector(hOld_[d],  (3 + d)     * blockSize);

        // neighbors
        if (nbrDofs)
        {
            mfem::Vector &eNbr = eOld_[d].FaceNbrData();
            mfem::Vector &hNbr = hOld_[d].FaceNbrData();
            MFEM_ASSERT(eNbr.Size() == nbrDofs, "e-nbr size mismatch");
            MFEM_ASSERT(hNbr.Size() == nbrDofs, "h-nbr size mismatch");

            inNew.SetVector(eNbr,      d      * blockSize + ndofs);
            inNew.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
        }
    }
#endif
    timerAssembleInNew.Stop();

    // ---------------------------
    // 3) Apply global operator A
    // ---------------------------
    timerApplyA.Start();
    out.SetSize(6 * ndofs);
    out.UseDevice(true);
    globalOperator_->Mult(inNew, out);
    timerApplyA.Stop();

    // ---------------------------
    // 4) Add forcing s(t_stage) (TFSF etc.), no source lifetime changes
    // ---------------------------
    timerTFSF.Start();
    if (auto *tfsf_space = srcmngr_.getGlobalTFSFSpace())
    {
        const int ndofs_tfsf = tfsf_space->GetNDofs();

        for (auto& source : srcmngr_.sources)
        {
            if (!source) { continue; }
            if (!dynamic_cast<TotalField*>(source.get())) { continue; }

#ifdef SEMBA_DGTD_ENABLE_CUDA
            const auto func_gpu = eval_time_var_field_gpu(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func_gpu);
#else
            auto func = evalTimeVarFunction(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
#endif
            mfem::Vector tempTFSF(assembledFunc.Size());
            tempTFSF.UseDevice(true);
            TFSFOperator_->Mult(assembledFunc, tempTFSF);

#ifdef SEMBA_DGTD_ENABLE_CUDA
            // Adds (-tempTFSF) into 'out'
            load_tfsf_into_out_vector_gpu(sub_to_parent_ids_, tempTFSF, out,
                                          ndofs, ndofs_tfsf);
#else
            for (int f : { E, H })
            {
                for (int d : { X, Y, Z })
                {
                    for (int v = 0; v < sub_to_parent_ids_.Size(); ++v)
                    {
                        const int outIdx  = (f * 3 + d) * ndofs + sub_to_parent_ids_[v];
                        const int tempIdx = (f * 3 + d) * ndofs_tfsf + v;
                        const double val = tempTFSF[tempIdx];
                        if (val != 0.0) { out[outIdx] -= val; }
                    }
                }
            }
#endif
        }
    }
    timerTFSF.Stop();

    timerTotal.Stop();

    // ---------------------------
    // 5) Periodic timing print
    // ---------------------------
    static double lastPrintTime = -1.0;
    const double printInterval = 0.1;  // in problem time units
    const double currentTime = GetTime();

    if (lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval)
    {
        std::cout << "Current time: " << currentTime << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " Mult total: "     << timerTotal.RealTime()        * 1000.0
                  << " ms, exchange: "   << timerExchange.RealTime()     * 1000.0
                  << " ms, assembleIn: " << timerAssembleInNew.RealTime()* 1000.0 << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " Mult applyA: "    << timerApplyA.RealTime()       * 1000.0
                  << " ms, tfsf: "       << timerTFSF.RealTime()         * 1000.0 << std::endl;

        lastPrintTime = (lastPrintTime < 0.0) ? currentTime : lastPrintTime + printInterval;
    }
}

void GlobalEvolution::ImplicitSolve(const double dt,
                                    const mfem::Vector& x,
                                    mfem::Vector& k)
{
    // Generic implicit stage: (I - dt * A) k = A * x + s(t_stage)
    // GetTime() is the stage time set by the integrator.
    mfem::StopWatch timerTotal, timerRHS, timerAx, timerSrc, timerSolve;

    timerTotal.Start();

    const int n = x.Size();
    const auto ndofs     = fes_.GetNDofs();
    const auto nbrDofs   = fes_.num_face_nbr_dofs;
    const auto blockSize = ndofs + nbrDofs;
    MFEM_ASSERT(n == 6 * ndofs, "ImplicitSolve: size mismatch");

    // ---- Helper: apply A (no sources), no printing here (used by GMRES) ----
    auto applyA_noPrint = [&](const mfem::Vector& u, mfem::Vector& Au)
    {
#ifdef SEMBA_DGTD_ENABLE_CUDA
        load_in_to_eh_gpu(u, eOld_, hOld_, ndofs);
#else
        for (int d = X; d <= Z; ++d)
        {
            for (int v = 0; v < ndofs; ++v)
            {
                eOld_[d][v] = u[d * ndofs + v];
                hOld_[d][v] = u[(3 + d) * ndofs + v];
            }
        }
#endif
        for (int d = X; d <= Z; ++d)
        {
            eOld_[d].ExchangeFaceNbrData();
            hOld_[d].ExchangeFaceNbrData();
        }

        mfem::Vector inNew(6 * blockSize);
        inNew.UseDevice(true);
#ifdef SEMBA_DGTD_ENABLE_CUDA
        load_eh_to_innew_gpu(u, inNew, ndofs, nbrDofs);
        load_nbr_to_innew_gpu(eOld_, hOld_, inNew, ndofs, nbrDofs);
#else
        for (int d = X; d <= Z; ++d)
        {
            inNew.SetVector(eOld_[d],       d      * blockSize);
            inNew.SetVector(hOld_[d],  (3 + d)     * blockSize);
            if (nbrDofs)
            {
                mfem::Vector &eNbr = eOld_[d].FaceNbrData();
                mfem::Vector &hNbr = hOld_[d].FaceNbrData();
                inNew.SetVector(eNbr,      d      * blockSize + ndofs);
                inNew.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
            }
        }
#endif
        Au.SetSize(6 * ndofs);
        Au.UseDevice(true);
        globalOperator_->Mult(inNew, Au);
    };

    // ---- Build RHS = A*x + s(t_stage), with timing breakdown ----
    timerRHS.Start();

    timerAx.Start();
    mfem::Vector rhs(n); rhs.UseDevice(true);
    applyA_noPrint(x, rhs);   // counts as A(x)
    timerAx.Stop();

    timerSrc.Start();
    {
        mfem::Vector s(n); s.UseDevice(true); s = 0.0;
        if (auto *tfsf_space = srcmngr_.getGlobalTFSFSpace())
        {
            const int ndofs_tfsf = tfsf_space->GetNDofs();

            for (auto& source : srcmngr_.sources)
            {
                if (!source) { continue; }
                if (!dynamic_cast<TotalField*>(source.get())) { continue; }

#ifdef SEMBA_DGTD_ENABLE_CUDA
                const auto func_gpu = eval_time_var_field_gpu(GetTime(), srcmngr_);
                mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func_gpu);
#else
                auto func = evalTimeVarFunction(GetTime(), srcmngr_);
                mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
#endif
                mfem::Vector tempTFSF(assembledFunc.Size());
                tempTFSF.UseDevice(true);
                TFSFOperator_->Mult(assembledFunc, tempTFSF);

#ifdef SEMBA_DGTD_ENABLE_CUDA
                // Adds (-tempTFSF) into 's'
                load_tfsf_into_out_vector_gpu(sub_to_parent_ids_, tempTFSF, s,
                                              ndofs, ndofs_tfsf);
#else
                for (int f : { E, H })
                {
                    for (int d : { X, Y, Z })
                    {
                        for (int v = 0; v < sub_to_parent_ids_.Size(); ++v)
                        {
                            const int outIdx  = (f * 3 + d) * ndofs + sub_to_parent_ids_[v];
                            const int tempIdx = (f * 3 + d) * ndofs_tfsf + v;
                            const double val = tempTFSF[tempIdx];
                            if (val != 0.0) { s[outIdx] -= val; }
                        }
                    }
                }
#endif
            }
        }
        rhs += s;
    }
    timerSrc.Stop();

    timerRHS.Stop();

    // ---- Define J(v) = v - dt * A v ----
    class JOp final : public mfem::Operator
    {
        std::function<void(const mfem::Vector&, mfem::Vector&)> applyA_;
        const double dt_;
        mutable mfem::Vector tmp_;
    public:
        JOp(int n,
            std::function<void(const mfem::Vector&, mfem::Vector&)> applyA,
            double dt)
        : mfem::Operator(n), applyA_(std::move(applyA)), dt_(dt), tmp_(n)
        { tmp_.UseDevice(true); }

        void Mult(const mfem::Vector& v, mfem::Vector& Jv) const override
        {
            applyA_(v, tmp_);   // tmp_ = A v
            Jv = v;             // Jv  = v
            Jv.Add(-dt_, tmp_); // Jv -= dt * A v
        }
    } J(n, std::function<void(const mfem::Vector&, mfem::Vector&)>(applyA_noPrint), dt);

    // ---- Solve ----
    mfem::GMRESSolver gmres;
    gmres.SetOperator(J);
    gmres.SetRelTol(1e-8);
    gmres.SetMaxIter(400);
    gmres.SetPrintLevel(0);
    // Optional: gmres.SetPreconditioner(P);

    timerSolve.Start();
    gmres.Mult(rhs, k);
    timerSolve.Stop();

    timerTotal.Stop();

    // ---------------------------
    // Periodic timing print
    // ---------------------------
    static double lastPrintTime = -1.0;
    const double printInterval = 0.1;  // in problem time units
    const double currentTime = GetTime();

    if (lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval)
    {
        std::cout << "Current time: " << currentTime << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " ImplicitSolve total: " << timerTotal.RealTime() * 1000.0
                  << " ms, rhs: " << timerRHS.RealTime() * 1000.0
                  << " ms A(x): " << timerAx.RealTime() * 1000.0
                  << " ms, source: " << timerSrc.RealTime() * 1000.0 << " ms)" << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " ImplicitSolve solve: " << timerSolve.RealTime() * 1000.0
                  << " ms, gmres iters: " << gmres.GetNumIterations() << std::endl;

        lastPrintTime = (lastPrintTime < 0.0) ? currentTime : lastPrintTime + printInterval;
    }
}

}
