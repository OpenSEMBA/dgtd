#include "GlobalEvolution.h"
#include "solver/SolverExtension.h"

#include <chrono>

namespace maxwell {

InteriorFaceConnectivityMaps getGlobalNodeID(const InteriorFaceConnectivityMaps& local_dof_ids, const GlobalConnectivity& global)
{
    InteriorFaceConnectivityMaps res;
    res.first.resize(local_dof_ids.first.size());
    res.second.resize(local_dof_ids.second.size());
    for (auto v{ 0 }; v < res.first.size(); v++) {
        res.first[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(local_dof_ids.first[v], local_dof_ids.second[v])));
        res.second[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(local_dof_ids.second[v], local_dof_ids.first[v])));
    }
    return res;
}

std::map<GeomTag, std::vector<NodePair>> GlobalEvolution::findSGBCDoFPairs()
{
    std::map<GeomTag, std::vector<NodePair>> res;
    
    auto attMap{ mapOriginalAttributes(*fes_.GetMesh()) };
    auto fec = dynamic_cast<const mfem::DG_FECollection*>(fes_.FEColl());
    GlobalConnectivity global = assembleGlobalConnectivityMap(*fes_.GetMesh(), fec);
    auto sgbc_marker = model_.getMarker(BdrCond::SGBC, true);
    if (sgbc_marker.Size() != 0){
        for (auto b = 0; b < model_.getMesh().GetNBE(); b++){
            if (sgbc_marker[model_.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
                const mfem::FaceElementTransformations* faceTrans;
                fes_.GetMesh()->FaceIsInterior(fes_.GetMesh()->GetFaceElementTransformations(fes_.GetMesh()->GetBdrElementFaceIndex(b))->ElementNo) ? faceTrans = fes_.GetMesh()->GetInternalBdrFaceTransformations(b) : faceTrans = fes_.GetMesh()->GetBdrFaceTransformations(b);
                if (fes_.GetMesh()->Dimension() == 1){
                    res[model_.getConstMesh().GetBdrAttribute(b)].emplace_back(faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder(), faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder() + 1);
                } else {
                    auto twoElemSubMesh{ assembleInteriorFaceSubMesh(*fes_.GetMesh(), *faceTrans, attMap) };
                    mfem::FiniteElementSpace subFES(&twoElemSubMesh, fec);
                    auto node_pair_global{ getGlobalNodeID(buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES), global)};
                    for (auto p = 0; p < node_pair_global.first.size(); p++){
                        res[model_.getConstMesh().GetBdrAttribute(b)].emplace_back(node_pair_global.first[p], node_pair_global.second[p]);
                    }
                }
            }
        }
    }
    return res;
}

FieldGridFuncs initSGBCHelperFields(const int size)
{
    FieldGridFuncs res;
    for (auto f : {E, H}){
        for (auto d : {X, Y , Z}){
            res.at(f).at(d).UseDevice(true);
            res.at(f).at(d).SetSize(size);
            res.at(f).at(d) = 0.0;
        }
    }
    return res;
}


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

    sgbc_pairs_ = findSGBCDoFPairs();
    for (const auto& [tag, pair_vector] : sgbc_pairs_){
		const auto& sbcps = model_.getSGBCProperties();
        for (auto p = 0; p < sbcps.size(); p++){
            for (auto t = 0; t < sbcps[p].geom_tags.size(); t++){
                if (tag == sbcps[p].geom_tags[t]){
                    auto wrap = SGBCWrapper::buildSGBCWrapper(sbcps[p]);
                    sgbcWrappers_.push_back(std::move(wrap));
                }
            }
        }
    }

    globalOperator_ = std::make_unique<mfem::SparseMatrix>(
        numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs(), 
        numberOfFieldComponents * numberOfMaxDimensions * (fes_.GetNDofs() + fes_.num_face_nbr_dofs)
    );

#ifdef SHOW_TIMER_INFORMATION
    if (Mpi::WorldRank() == 0){
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }
#endif


    Probes probes; // Local probes instance for operator construction
    if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {

        srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));

        auto globalTFSFfes = srcmngr_.getGlobalTFSFSpace();
        auto tfsfMesh = globalTFSFfes->GetMesh();

        Model tfsfModel = Model(*tfsfMesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
        
        ProblemDescription tfsfpd(tfsfModel, probes, srcmngr_.sources, opts_);
        DGOperatorFactory<FiniteElementSpace> tfsfops(tfsfpd, *globalTFSFfes);

        TFSFOperator_ = tfsfops.buildTFSFGlobalOperator();

        auto src_sm = static_cast<mfem::SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
        mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), tfsf_sub_to_parent_ids_);
    }

    if (model_.getSGBCToMarker().find(BdrCond::SGBC) != model_.getSGBCToMarker().end()) {

        srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getSGBCToMarker().at(BdrCond::SGBC));

        auto globalSGBCfes = srcmngr_.getGlobalTFSFSpace();
        auto sgbcMesh = globalSGBCfes->GetMesh();

        Model sgbcModel = Model(*sgbcMesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
        
        ProblemDescription sgbcpd(sgbcModel, probes, srcmngr_.sources, opts_);
        DGOperatorFactory<FiniteElementSpace> sgbcops(sgbcpd, *globalSGBCfes);

        SGBCOperator_ = sgbcops.buildTFSFGlobalOperator();

        auto src_sm = static_cast<mfem::SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
        mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), sgbc_sub_to_parent_ids_);
    }

    ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
    DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, fes_);

    globalOperator_ = dgops.buildGlobalOperator();

}

const mfem::Vector buildSingleVectorTFSFFunc(const FieldGridFuncs& func)
{
    mfem::Vector res(6 * func[E][X].Size());
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

void assertVectorOnDevice(const mfem::Vector &v, const std::string &name)
{
    mfem::MemoryType mem_type = v.GetMemory().GetMemoryType();
    if (mem_type != mfem::MemoryType::DEVICE)
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
    mfem::StopWatch timerTotal, timerExchange, timerAssembleInNew, timerApplyA, timerTFSF, timerSGBC;

    auto current_debug_time = GetTime();

    timerTotal.Start();

    const auto ndofs     = fes_.GetNDofs();
    const auto nbrDofs   = fes_.num_face_nbr_dofs;
    const auto blockSize = ndofs + nbrDofs;

    // 1) Load locals and exchange neighbors
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

    // 2) Assemble packed input [locals|neighbors]
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

    // 3) Apply global operator A
    timerApplyA.Start();
    out.SetSize(6 * ndofs);
    out.UseDevice(true);
    globalOperator_->Mult(inNew, out);
    timerApplyA.Stop();

    // 4) Add forcing s(t_stage) (TFSF etc.), gated by final time
    timerTFSF.Start();
    if (auto *tfsf_space = srcmngr_.getGlobalTFSFSpace())
    {
        const double t_stage = GetTime();
        const double t_off   = opts_.tfsf_cutoff_time;

        // Gate: if tfsfFinalTime is set and we are at/after it, skip TFSF assembly entirely.
        const bool tfsf_active = (t_off == 0.0) || (t_stage < t_off);

        if (tfsf_active){
            const int ndofs_tfsf = tfsf_space->GetNDofs();

            for (auto& source : srcmngr_.sources){
                if (!source) { 
                    continue; 
                }
                if (!dynamic_cast<TotalField*>(source.get())) { 
                    continue; 
                }

#ifdef SEMBA_DGTD_ENABLE_CUDA
                const auto func_gpu = eval_time_var_field_gpu(t_stage, srcmngr_);
                mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func_gpu);
#else
                auto func = evalTimeVarFunction(t_stage, srcmngr_);
                mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
#endif
                mfem::Vector tempTFSF(assembledFunc.Size());
                tempTFSF.UseDevice(true);
                TFSFOperator_->Mult(assembledFunc, tempTFSF);

#ifdef SEMBA_DGTD_ENABLE_CUDA
                // Adds (-tempTFSF) into 'out'
                load_tfsf_into_out_vector_gpu(tfsf_sub_to_parent_ids_, tempTFSF, out,
                                              ndofs, ndofs_tfsf);
#else
                for (int f : { E, H }){
                    for (int d : { X, Y, Z }){
                        for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v){
                            const int outIdx  = (f * 3 + d) * ndofs + tfsf_sub_to_parent_ids_[v];
                            const int tempIdx = (f * 3 + d) * ndofs_tfsf + v;
                            const double val  = tempTFSF[tempIdx];
                            out[outIdx] -= val;
                        }
                    }
                }
#endif
            }
        }
    }
    timerTFSF.Stop();

    timerSGBC.Start();
    for (const auto& [tag, pairs] : sgbc_pairs_){
        for (auto w = 0; w < sgbcWrappers_.size(); w++){
            if (tag == sgbcWrappers_[w]->getProperties().geom_tags.at(0)){
                for (auto p = 0; p < pairs.size(); p++){
                    
                    sgbcWrappers_[w]->updateFieldsWithGlobal(eOld_, hOld_, pairs[p]);
                    sgbcWrappers_[w]->solve(GetTime(), GetTime() - sgbcWrappers_[w]->getOldTime());
                    sgbcWrappers_[w]->setOldTime(GetTime());

                    const auto sgbc_vec_size = SGBCOperator_->Height() / 6;
                    auto sgbc_fields = initSGBCHelperFields(sgbc_vec_size);
                    sgbcWrappers_[w]->getSGBCFields(sgbc_sub_to_parent_ids_, pairs[p], sgbc_fields);
                    auto sgbc_vector = buildSingleVectorTFSFFunc(sgbc_fields);
                    mfem::Vector tempSGBC(sgbc_vector.Size());
                    tempSGBC.UseDevice(true);

                    SGBCOperator_->Mult(sgbc_vector, tempSGBC);

                    auto dir_size = tempSGBC.Size() / 6;
                    for (int f : { E, H }){
                        for (int d : { X, Y, Z }){
                            for (auto v = 0; v < dir_size / 2; v++){
                                tempSGBC[(f * 3 + d) * dir_size + dir_size - 1 - v] = tempSGBC[(f * 3 + d) * dir_size + v];
                                // tempSGBC[(f * 3 + d) * dir_size + v] -= tempSGBC[(f * 3 + d) * dir_size + dir_size - 1 - v];
                            }
                        }
                    }

                    for (int f : { E, H }){
                        for (int d : { X, Y, Z }){
                            for (int v = 0; v < sgbc_sub_to_parent_ids_.Size(); ++v){
                                const int out_idx = (f * 3 + d) * ndofs + sgbc_sub_to_parent_ids_[v];
                                const int temp_idx = (f * 3 + d) * sgbc_vec_size + v;
                                if (v < sgbc_sub_to_parent_ids_.Size() / 2){
                                    out[out_idx] += tempSGBC[temp_idx];
                                } else {
                                    out[out_idx] += tempSGBC[temp_idx];
                                }
                            }
                        }
                    }

                }
            }
        }
    }
    timerSGBC.Stop();

    timerTotal.Stop();

    // 5) Periodic timing print
    static double lastPrintTime = -1.0;
    const double printInterval = 0.1;  // in problem time units
    const double currentTime = GetTime();

    if (lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval)
    {
        std::cout << "Current time: " << currentTime << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " Mult total: "     << timerTotal.RealTime()         * 1000.0
                  << " ms, exchange: "   << timerExchange.RealTime()      * 1000.0
                  << " ms, assembleIn: " << timerAssembleInNew.RealTime() * 1000.0 << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " Mult applyA: "    << timerApplyA.RealTime()        * 1000.0
                  << " ms, tfsf: "       << timerTFSF.RealTime()          * 1000.0 
                  << " ms, sgbc: "       << timerSGBC.RealTime()          * 1000.0 << std::endl;

        lastPrintTime = (lastPrintTime < 0.0) ? currentTime : lastPrintTime + printInterval;
    }
}

void GlobalEvolution::ImplicitSolve(const double dt,
                                    const mfem::Vector& x,
                                    mfem::Vector& k)
{
    // (I - dt * A) k = A * x + s(t_stage), with GetTime() already set to stage time.
    mfem::StopWatch timerTotal, timerRHS, timerAx, timerSrc, timerSolve;
    timerTotal.Start();

    const int n         = x.Size();
    const int ndofs     = fes_.GetNDofs();
    const int nbrDofs   = fes_.num_face_nbr_dofs;
    const int blockSize = ndofs + nbrDofs;
    MFEM_ASSERT(n == 6 * ndofs, "ImplicitSolve: size mismatch");

    // Reusable work buffers: eliminate per-matvec allocations/copies inside GMRES.
    struct ApplyAWork {
        mfem::Vector inNew;  // packed [locals|neighbors]
        mfem::Vector Au;     // output buffer
        ApplyAWork(int packed, int out) : inNew(packed), Au(out) {
            inNew.UseDevice(true);
            Au.UseDevice(true);
        }
    } work(6 * blockSize, 6 * ndofs);

    // Matrix-free apply of A (no sources), no printing inside (used by GMRES).
    auto applyA_noPrint = [&](const mfem::Vector& u, mfem::Vector& Au)
    {
#ifdef SEMBA_DGTD_ENABLE_CUDA
        load_in_to_eh_gpu(u, eOld_, hOld_, ndofs);
#else
        for (int d = X; d <= Z; ++d) {
            for (int v = 0; v < ndofs; ++v) {
                eOld_[d][v] = u[d * ndofs + v];
                hOld_[d][v] = u[(3 + d) * ndofs + v];
            }
        }
#endif
        for (int d = X; d <= Z; ++d) {
            eOld_[d].ExchangeFaceNbrData();
            hOld_[d].ExchangeFaceNbrData();
        }

#ifdef SEMBA_DGTD_ENABLE_CUDA
        load_eh_to_innew_gpu(u, work.inNew, ndofs, nbrDofs);
        load_nbr_to_innew_gpu(eOld_, hOld_, work.inNew, ndofs, nbrDofs);
#else
        for (int d = X; d <= Z; ++d) {
            work.inNew.SetVector(eOld_[d],       d      * blockSize);
            work.inNew.SetVector(hOld_[d],  (3 + d)     * blockSize);
            if (nbrDofs) {
                mfem::Vector &eNbr = eOld_[d].FaceNbrData();
                mfem::Vector &hNbr = hOld_[d].FaceNbrData();
                work.inNew.SetVector(eNbr,      d      * blockSize + ndofs);
                work.inNew.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
            }
        }
#endif

        globalOperator_->Mult(work.inNew, work.Au);
        Au = work.Au; // alias-safe copy
    };

    // Build RHS = A*x + s(t_stage) with timing breakdown.
    timerRHS.Start();

    timerAx.Start();
    mfem::Vector rhs(n); rhs.UseDevice(true);
    applyA_noPrint(x, rhs);
    timerAx.Stop();

    timerSrc.Start();
    {
        mfem::Vector s(n); s.UseDevice(true); s = 0.0;

        if (auto *tfsf_space = srcmngr_.getGlobalTFSFSpace()) {
            const double t_stage = GetTime();
            const double t_off   = opts_.tfsf_cutoff_time;

            // Gate: if tfsfFinalTime is set and we are at/after it, skip TFSF assembly entirely.
            const bool tfsf_active = (t_off == 0.0) || (t_stage < t_off);

            if (tfsf_active)
            {
                const int ndofs_tfsf = tfsf_space->GetNDofs();

                for (auto& source : srcmngr_.sources) {
                    if (!source) { continue; }
                    if (!dynamic_cast<TotalField*>(source.get())) { continue; }

#ifdef SEMBA_DGTD_ENABLE_CUDA
                    const auto func_gpu = eval_time_var_field_gpu(t_stage, srcmngr_);
                    mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func_gpu);
#else
                    auto func = evalTimeVarFunction(t_stage, srcmngr_);
                    mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
#endif
                    mfem::Vector tempTFSF(assembledFunc.Size());
                    tempTFSF.UseDevice(true);
                    TFSFOperator_->Mult(assembledFunc, tempTFSF);

#ifdef SEMBA_DGTD_ENABLE_CUDA
                    // Adds (-tempTFSF) into 's'
                    load_tfsf_into_out_vector_gpu(tfsf_sub_to_parent_ids_, tempTFSF, s,
                                                  ndofs, ndofs_tfsf);
#else
                    for (int f : { E, H }) {
                        for (int d : { X, Y, Z }) {
                            for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v) {
                                const int outIdx  = (f * 3 + d) * ndofs + tfsf_sub_to_parent_ids_[v];
                                const int tempIdx = (f * 3 + d) * ndofs_tfsf + v;
                                const double val = tempTFSF[tempIdx];
                                if (val != 0.0) { s[outIdx] -= val; }
                            }
                        }
                    }
#endif
                }
            }
        }

        rhs += s;
    }
    timerSrc.Stop();

    timerRHS.Stop();

    // Linear operator: J(v) = v - dt * A v
    class JOp final : public mfem::Operator {
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
            applyA_(v, tmp_);   // tmp = A v
            Jv = v;             // Jv  = v
            Jv.Add(-dt_, tmp_); // Jv -= dt * A v
        }
    } J(n, std::function<void(const mfem::Vector&, mfem::Vector&)>(applyA_noPrint), dt);

    // GMRES solve without preconditioner.
    mfem::GMRESSolver gmres;
    gmres.SetOperator(J);
    gmres.SetRelTol(1e-8);
    gmres.SetMaxIter(400);
    gmres.SetKDim(60);
    gmres.SetPrintLevel(0);

    // Initial guess: reuse current k if sized; else use rhs.
    if (k.Size() == n) {
        // keep caller-provided guess
    } else {
        k.SetSize(n); k.UseDevice(true); k = rhs;
    }

    timerSolve.Start();
    gmres.Mult(rhs, k);
    timerSolve.Stop();

    timerTotal.Stop();

    // Periodic timing print (problem time).
    static double lastPrintTime = -1.0;
    const double printInterval = 0.1;
    const double currentTime   = GetTime();

    if (lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval) {
        std::cout << "Current time: " << currentTime << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " ImplicitSolve total: " << timerTotal.RealTime() * 1000.0
                  << " ms, rhs: " << timerRHS.RealTime() * 1000.0
                  << " ms (A(x): " << timerAx.RealTime() * 1000.0
                  << " ms, source: " << timerSrc.RealTime() * 1000.0 << " ms)" << std::endl;
        std::cout << "Rank " << Mpi::WorldRank()
                  << " ImplicitSolve solve: " << timerSolve.RealTime() * 1000.0
                  << " ms, gmres iters: " << gmres.GetNumIterations() << std::endl;

        lastPrintTime = (lastPrintTime < 0.0) ? currentTime : lastPrintTime + printInterval;
    }
}

}