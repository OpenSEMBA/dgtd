#include "GlobalEvolution.h"
#include "solver/SolverExtension.h"

#include <chrono>
#include <unordered_set>

namespace maxwell {

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

    std::map<GeomTag, std::vector<NodePair>> raw_pairs;
    {
        auto attMap{ mapOriginalAttributes(*fes_.GetMesh()) };
        auto fec = dynamic_cast<const mfem::DG_FECollection*>(fes_.FEColl());
        
        auto sgbc_marker = model_.getMarker(BdrCond::SGBC, true);
        if (sgbc_marker.Size() != 0){
            for (auto b = 0; b < model_.getMesh().GetNBE(); b++){
                if (sgbc_marker[model_.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
                    const mfem::FaceElementTransformations* faceTrans;
                    fes_.GetMesh()->FaceIsInterior(fes_.GetMesh()->GetFaceElementTransformations(fes_.GetMesh()->GetBdrElementFaceIndex(b))->ElementNo) ? faceTrans = fes_.GetMesh()->GetInternalBdrFaceTransformations(b) : faceTrans = fes_.GetMesh()->GetBdrFaceTransformations(b);
                    
                    if (fes_.GetMesh()->Dimension() == 1){
                        raw_pairs[model_.getConstMesh().GetBdrAttribute(b)].emplace_back(
                            faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder(), 
                            faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder() + 1
                        );
                    } else {
                        auto twoElemSubMesh{ assembleInteriorFaceSubMesh(*fes_.GetMesh(), *faceTrans, attMap) };
                        mfem::FiniteElementSpace subFES(&twoElemSubMesh, fec);
                        
                        auto node_pair_local = buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES);
                        
                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            raw_pairs[model_.getConstMesh().GetBdrAttribute(b)].emplace_back(
                                node_pair_local.first[p], 
                                node_pair_local.second[p]
                            );
                        }
                    }
                }
            }
        }
    }

    mfem::Array<int> temp_dof_to_elem(fes_.GetNDofs());
    temp_dof_to_elem = -1;
    {
        const mfem::Table& elem_dof = fes_.GetElementToDofTable();
        for (int el = 0; el < elem_dof.Size(); el++) {
            const int* dofs = elem_dof.GetRow(el);
            const int n_dofs = elem_dof.RowSize(el);
            for (int i = 0; i < n_dofs; i++) {
                temp_dof_to_elem[dofs[i]] = el;
            }
        }
    }

    for (const auto& [tag, pairs] : raw_pairs){
        const auto& sbcps = model_.getSGBCProperties();
        for (auto p = 0; p < sbcps.size(); p++){
            for (auto t = 0; t < sbcps[p].geom_tags.size(); t++){
                if (tag == sbcps[p].geom_tags[t]){
                    
                    auto wrap = SGBCWrapper::buildSGBCWrapper(sbcps[p]);
                    
                    int state_size = wrap->getStateSize();
                    for(const auto& pair : pairs) {
                        SGBCState s;
                        s.global_pair = pair;
                        
                        s.element_pair.first = temp_dof_to_elem[pair.first];
                        s.element_pair.second = temp_dof_to_elem[pair.second];
                        
                        if (s.element_pair.first == -1 || s.element_pair.second == -1) {
                            std::cerr << "Error: Could not find Element ID for SGBC Node Pair: " 
                                      << pair.first << ", " << pair.second << std::endl;
                        }

                        s.init(state_size);
                        sgbc_states_[tag].push_back(s);
                    }
                    
                    sgbcWrappers_.push_back(std::move(wrap));
                }
            }
        }
    }

    globalOperator_ = std::make_unique<mfem::SparseMatrix>(
        numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs(), 
        numberOfFieldComponents * numberOfMaxDimensions * (fes_.GetNDofs() + fes_.num_face_nbr_dofs)
    );

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
        mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), tfsf_sub_to_parent_ids_);
    }

    if (model_.getSGBCToMarker().find(BdrCond::SGBC) != model_.getSGBCToMarker().end()) {
        auto sgbc_sm = TotalFieldScatteredFieldSubMesher(model_.getConstMesh(), model_.getSGBCToMarker().at(BdrCond::SGBC));
        auto global_sm_fes = std::make_unique<FiniteElementSpace>(sgbc_sm.getGlobalTFSFSubMesh(), fes_.FEColl());
        SGBCndofs_ = global_sm_fes->GetNDofs();
        auto sgbc_mesh = global_sm_fes->GetMesh();
        auto src_sm = static_cast<mfem::SubMesh*>(sgbc_mesh);
        mfem::SubMeshUtils::BuildVdofToVdofMap(*global_sm_fes, fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), sgbc_sub_to_parent_ids_);
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
    if (mem_type != mfem::MemoryType::DEVICE) {
        mfem::out << "Warning: Vector '" << name << "' latest data is NOT on device!\n";
    } else {
        mfem::out << "Vector '" << name << "' latest data is on DEVICE.\n";
    }
}

void GlobalEvolution::advanceSGBCs(double time, double dt, 
                                   const std::array<mfem::ParGridFunction, 3>& e, 
                                   const std::array<mfem::ParGridFunction, 3>& h)
{
    for (auto& [tag, states] : sgbc_states_){
        for (auto& w : sgbcWrappers_){
            if (tag == w->getProperties().geom_tags.at(0)){
                
                for (auto& state : states){
                    w->loadState(state);
                    w->updateFieldsWithGlobal(e, h, state);
                    w->solve(time, dt);
                    w->saveState(state);
                }
                w->setOldTime(time); 
            }
        }
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

    timerExchange.Start();
#ifdef SEMBA_DGTD_ENABLE_CUDA
    load_in_to_eh_gpu(in, eOld_, hOld_, ndofs);
#else
    for (int d = X; d <= Z; ++d) {
        for (int v = 0; v < ndofs; ++v) {
            eOld_[d][v] = in[d * ndofs + v];
            hOld_[d][v] = in[(3 + d) * ndofs + v];
        }
    }
#endif
    for (int d = X; d <= Z; ++d) {
        eOld_[d].ExchangeFaceNbrData();
        hOld_[d].ExchangeFaceNbrData();
    }
    timerExchange.Stop();

    timerAssembleInNew.Start();
    mfem::Vector inNew(6 * blockSize);
    inNew.UseDevice(true);

#ifdef SEMBA_DGTD_ENABLE_CUDA
    load_eh_to_innew_gpu(in, inNew, ndofs, nbrDofs);
    load_nbr_to_innew_gpu(eOld_, hOld_, inNew, ndofs, nbrDofs);
#else
    for (int d = X; d <= Z; ++d) {
        inNew.SetVector(eOld_[d],       d      * blockSize);
        inNew.SetVector(hOld_[d],  (3 + d)     * blockSize);
        if (nbrDofs) {
            mfem::Vector &eNbr = eOld_[d].FaceNbrData();
            mfem::Vector &hNbr = hOld_[d].FaceNbrData();
            inNew.SetVector(eNbr,      d      * blockSize + ndofs);
            inNew.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
        }
    }
#endif
    timerAssembleInNew.Stop();

    timerApplyA.Start();
    out.SetSize(6 * ndofs);
    out.UseDevice(true);
    globalOperator_->Mult(inNew, out);
    timerApplyA.Stop();

    timerTFSF.Start();
    if (auto *tfsf_space = srcmngr_.getGlobalTFSFSpace())
    {
        const double t_stage = GetTime();
        const double t_off   = opts_.tfsf_cutoff_time;
        const bool tfsf_active = (t_off == 0.0) || (t_stage < t_off);

        if (tfsf_active){
            const int ndofs_tfsf = tfsf_space->GetNDofs();
            for (auto& source : srcmngr_.sources){
                if (!source) continue; 
                if (!dynamic_cast<TotalField*>(source.get())) continue; 

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
                load_tfsf_into_out_vector_gpu(tfsf_sub_to_parent_ids_, tempTFSF, out, ndofs, ndofs_tfsf);
#else
                for (int f : { E, H }){
                    for (int d : { X, Y, Z }){
                        for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v){
                            const int outIdx  = (f * 3 + d) * ndofs + tfsf_sub_to_parent_ids_[v];
                            const int tempIdx = (f * 3 + d) * ndofs_tfsf + v;
                            out[outIdx] -= tempTFSF[tempIdx];
                        }
                    }
                }
#endif
            }
        }
    }
    timerTFSF.Stop();

    // 5) SGBC Coupling via Flux Injection
    timerSGBC.Start();
    for (const auto& [tag, states] : sgbc_states_){
        for (auto& w : sgbcWrappers_){
            if (tag == w->getProperties().geom_tags.at(0)){
                
                auto sgbc_fields = initSGBCHelperFields(SGBCndofs_);
                
                for (const auto& state : states) {
                    
                    w->getSGBCFields(sgbc_sub_to_parent_ids_, state, sgbc_fields);
                    
                    const int global_node_left = state.global_pair.first;  
                    const int global_node_right = state.global_pair.second;
                    
                    const int sub_idx_L = sgbc_sub_to_parent_ids_.Find(global_node_left);
                    const int sub_idx_R = sgbc_sub_to_parent_ids_.Find(global_node_right);

                    if (sub_idx_L == -1 || sub_idx_R == -1) continue; 

                    // --- Mass Matrix Scaling Logic ---
                    // [UPDATED] Read Element IDs directly from the SGBCState
                    int el_idx_L = state.element_pair.first;
                    int el_idx_R = state.element_pair.second;

                    // Robust: Use element order (assumes pairs share order, or look up per element)
                    int order_L = fes_.GetElementOrder(el_idx_L); 
                    int order_R = fes_.GetElementOrder(el_idx_R); 

                    // Get exact volumes from the mesh using the stored element IDs
                    double vol_L = fes_.GetMesh()->GetElementVolume(el_idx_L);
                    double vol_R = fes_.GetMesh()->GetElementVolume(el_idx_R);

                    double scale_L = (double)(order_L * (order_L + 1)) / vol_L;
                    double scale_R = (double)(order_R * (order_R + 1)) / vol_R;

                    for (auto d : {X, Y, Z}) {
                        // --- Left Interface ---
                        double Eg_L = eOld_[d][global_node_left];
                        double Hg_L = hOld_[d][global_node_left];
                        double Es_L = sgbc_fields[E][d][sub_idx_L]; 
                        double Hs_L = sgbc_fields[H][d][sub_idx_L]; 

                        double flux_E_L = ((Hs_L - Hg_L) + (Es_L - Eg_L)) * 0.5;
                        double flux_H_L = ((Es_L - Eg_L) + (Hs_L - Hg_L)) * 0.5;

                        out[(E * 3 + d) * ndofs + global_node_left] += (flux_E_L * scale_L);
                        out[(H * 3 + d) * ndofs + global_node_left] += (flux_H_L * scale_L);

                        // --- Right Interface ---
                        double Eg_R = eOld_[d][global_node_right];
                        double Hg_R = hOld_[d][global_node_right];
                        double Es_R = sgbc_fields[E][d][sub_idx_R]; 
                        double Hs_R = sgbc_fields[H][d][sub_idx_R]; 

                        double flux_E_R = -((Hs_R - Hg_R) - (Es_R - Eg_R)) * 0.5;
                        double flux_H_R = -((Es_R - Eg_R) - (Hs_R - Hg_R)) * 0.5;

                        out[(E * 3 + d) * ndofs + global_node_right] += (flux_E_R * scale_R);
                        out[(H * 3 + d) * ndofs + global_node_right] += (flux_H_R * scale_R);
                    }
                }
            }
        }
    }
    timerSGBC.Stop();

    timerTotal.Stop();

    static double lastPrintTime = -1.0;
    const double printInterval = 0.1; 
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
    mfem::StopWatch timerTotal, timerRHS, timerAx, timerSrc, timerSolve;
    timerTotal.Start();

    const int n         = x.Size();
    const int ndofs     = fes_.GetNDofs();
    const int nbrDofs   = fes_.num_face_nbr_dofs;
    const int blockSize = ndofs + nbrDofs;
    MFEM_ASSERT(n == 6 * ndofs, "ImplicitSolve: size mismatch");

    struct ApplyAWork {
        mfem::Vector inNew;  
        mfem::Vector Au;     
        ApplyAWork(int packed, int out) : inNew(packed), Au(out) {
            inNew.UseDevice(true);
            Au.UseDevice(true);
        }
    } work(6 * blockSize, 6 * ndofs);

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
        Au = work.Au; 
    };

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
                    load_tfsf_into_out_vector_gpu(tfsf_sub_to_parent_ids_, tempTFSF, s, ndofs, ndofs_tfsf);
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
            applyA_(v, tmp_);   
            Jv = v;             
            Jv.Add(-dt_, tmp_); 
        }
    } J(n, std::function<void(const mfem::Vector&, mfem::Vector&)>(applyA_noPrint), dt);

    mfem::GMRESSolver gmres;
    gmres.SetOperator(J);
    gmres.SetRelTol(1e-8);
    gmres.SetMaxIter(400);
    gmres.SetKDim(60);
    gmres.SetPrintLevel(0);

    if (k.Size() == n) {
    } else {
        k.SetSize(n); k.UseDevice(true); k = rhs;
    }

    timerSolve.Start();
    gmres.Mult(rhs, k);
    timerSolve.Stop();

    timerTotal.Stop();

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