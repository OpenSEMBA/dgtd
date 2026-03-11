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
        auto mesh = fes_.GetMesh();
        auto const_mesh = &model_.getConstMesh();

        auto sgbc_marker_interior = model_.getMarker(BdrCond::SGBC, true);
        auto sgbc_marker_boundary = model_.getMarker(BdrCond::SGBC, false);

        if (sgbc_marker_interior.Size() != 0 || sgbc_marker_boundary.Size() != 0){
            for (auto b = 0; b < model_.getMesh().GetNBE(); b++){
                int bdr_attr = const_mesh->GetBdrAttribute(b);
                bool is_interior = false, is_boundary = false;

                if (sgbc_marker_interior.Size() != 0 && sgbc_marker_interior[bdr_attr - 1] == 1) {
                    is_interior = true;
                }
                if (sgbc_marker_boundary.Size() != 0 && sgbc_marker_boundary[bdr_attr - 1] == 1) {
                    is_boundary = true;
                }

                if (!is_interior && !is_boundary) continue;

                const mfem::FaceElementTransformations* faceTrans;
                mesh->FaceIsInterior(mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(b))->ElementNo) ?
                    faceTrans = mesh->GetInternalBdrFaceTransformations(b) :
                    faceTrans = mesh->GetBdrFaceTransformations(b);

                if (is_interior) {
                    if (mesh->Dimension() == 1){
                        raw_pairs[bdr_attr].emplace_back(
                            faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder(),
                            faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder() + 1
                        );
                    } else {
                        auto twoElemSubMesh{ assembleInteriorFaceSubMesh(*mesh, *faceTrans, attMap) };
                        mfem::FiniteElementSpace subFES(&twoElemSubMesh, fec);
                        auto node_pair_local = buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES);

                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            raw_pairs[bdr_attr].emplace_back(
                                node_pair_local.first[p],
                                node_pair_local.second[p]
                            );
                        }
                    }
                }

                if (is_boundary) {
                    if (mesh->Dimension() == 1){
                        if(faceTrans->Elem1No == 0) {
                            raw_pairs[bdr_attr].emplace_back(0, -1);
                        }
                        else {
                            raw_pairs[bdr_attr].emplace_back(faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder(), -1);
                        }
                    } else {
                        auto elementSubMesh{ assembleBoundaryFaceSubMesh(*mesh, *faceTrans, attMap)};
                        mfem::FiniteElementSpace subFES(&elementSubMesh, fec);
                        auto node_pair_local = buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES);
                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            raw_pairs[bdr_attr].emplace_back(
                                node_pair_local.first[p],
                                -1
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

    MPI_Barrier(MPI_COMM_WORLD);

    for (const auto& [tag, pairs] : raw_pairs){
        const auto& sbcps = model_.getSGBCProperties();
        for (auto p = 0; p < sbcps.size(); p++){
            for (auto t = 0; t < sbcps[p].geom_tags.size(); t++){
                if (tag == sbcps[p].geom_tags[t]){

                    auto wrap = SGBCWrapper::buildSGBCWrapper(sbcps[p]);
                    SGBCWrapper* wrap_ptr = wrap.get();

                    int state_size = wrap->getStateSize();
                    for(const auto& pair : pairs) {
                        SGBCState s;
                        s.global_pair = pair;

                        s.element_pair.first = temp_dof_to_elem[pair.first];
                        if (pair.second != -1){
                            s.element_pair.second = temp_dof_to_elem[pair.second];
                        }
                        else {
                            s.element_pair.second = -1;
                        }

                        if (s.element_pair.first == -1) {
                            std::cerr << "Error: Could not find Element ID for SGBC Node Pair: "
                                      << pair.first << ", " << pair.second << std::endl;
                        }

                        s.init(state_size);
                        sgbc_states_[tag].push_back(s);
                    }

                    sgbcWrappers_.push_back(std::move(wrap));
                    sgbc_wrapper_map_[tag] = wrap_ptr;
                }
            }
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

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
        std::cout << "TFSF elements in rank " << std::to_string(Mpi::WorldRank()) << ": " << globalTFSFfes->GetNE() << std::endl;
        if (globalTFSFfes->GetNE() <= 0){
            std::cout << "Warning! No TFSF elements on rank but it is defined to have, simulation may fail. Use less processes." << std::endl;
        }
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

        // --- NEW PREVENTATIVE FIX ---
        // Initialize all submesh boundary attributes to a safe dummy value (>0) 
        // to prevent `attr - 1 = -2` crash on the artificial cut faces.
        for (int c_be = 0; c_be < src_sm->GetNBE(); c_be++) {
            src_sm->SetBdrAttribute(c_be, 1); 
        }
        // ----------------------------

        auto parent_mesh = fes_.GetMesh();
        auto parent_f2bdr_map = parent_mesh->GetFaceToBdrElMap();
        auto child_f2bdr_map = src_sm->GetFaceToBdrElMap();
        auto face_map = mfem::SubMeshUtils::BuildFaceMap(*parent_mesh, *src_sm, src_sm->GetParentElementIDMap());

        GeomTagToBoundary sgbc_bdr;
        const auto& sgbc_marker = model_.getSGBCToMarker().at(BdrCond::SGBC);

        // Locate the pure boundary elements and transfer ONLY the SGBC tags to the submesh
        for (int be = 0; be < parent_mesh->GetNBE(); be++) {
            int parent_tag = parent_mesh->GetBdrAttribute(be);
            if (sgbc_marker[parent_tag - 1] == 1) {
                int p_face = parent_f2bdr_map.Find(be);
                if (p_face != -1) {
                    int c_face = face_map.Find(p_face);
                    if (c_face != -1) {
                        int c_be = child_f2bdr_map[c_face];
                        if (c_be != -1) {
                            src_sm->SetBdrAttribute(c_be, parent_tag);
                            sgbc_bdr[parent_tag] = BdrCond::SGBC;
                        }
                    }
                }
            }
        }
        
        src_sm->SetAttributes(); // Force MFEM to rebuild the bdr_attributes array

        Model sgbc_model = Model(*sgbc_mesh, GeomTagToMaterialInfo(), 
                                 GeomTagToBoundaryInfo(sgbc_bdr, GeomTagToInteriorBoundary{}));
        ProblemDescription sgbc_pd(sgbc_model, probes, srcmngr_.sources, opts_);
        DGOperatorFactory<FiniteElementSpace> sgbc_ops(sgbc_pd, *global_sm_fes);
        SGBCOperator_ = sgbc_ops.buildSGBCGlobalOperator();
        
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
    const int size = func[E][X].Size();
    for (auto f : { E, H }) {
        for (auto d : { X, Y, Z }) {
            for (auto v{ 0 }; v < size; v++) {
                res[v + (f * 3 + d) * size] = func[f][d][v];
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
        auto it = sgbc_wrapper_map_.find(tag);
        if (it == sgbc_wrapper_map_.end()) continue;

        auto w = it->second;
        for (auto& state : states){
            w->loadState(state);
            w->updateFieldsWithGlobal(e, h, state);
            w->solve(time, dt);
            w->saveState(state);
        }
        w->setOldTime(time);
    }
}

void GlobalEvolution::applyTFSFSourceToVector(double t_stage, int ndofs, int ndofs_tfsf,
                                              mfem::Vector& result_vector, bool check_zero) const
{
    if (auto *tfsf_space = srcmngr_.getGlobalTFSFSpace())
    {
        const double t_off   = opts_.tfsf_cutoff_time;
        const bool tfsf_active = (t_off == 0.0) || (t_stage < t_off);

        if (tfsf_active){
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
                load_tfsf_into_out_vector_gpu(tfsf_sub_to_parent_ids_, tempTFSF, result_vector, ndofs, ndofs_tfsf);
#else
                for (int f : { E, H }){
                    for (int d : { X, Y, Z }){
                        for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v){
                            const int outIdx  = (f * 3 + d) * ndofs + tfsf_sub_to_parent_ids_[v];
                            const int tempIdx = (f * 3 + d) * ndofs_tfsf + v;
                            const double val = tempTFSF[tempIdx];
                            if (!check_zero || val != 0.0) {
                                result_vector[outIdx] -= val;
                            }
                        }
                    }
                }
#endif
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
    }
    if (nbrDofs) {
        for (int d = X; d <= Z; ++d) {
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
    applyTFSFSourceToVector(GetTime(), ndofs, (srcmngr_.getGlobalTFSFSpace() ? srcmngr_.getGlobalTFSFSpace()->GetNDofs() : 0), out, false);
    timerTFSF.Stop();

    // 5) SGBC Coupling via Flux Injection
    timerSGBC.Start();

    if (sgbcWrappers_.size() != 0){
        // Reuse or initialize SGBC helper fields
        if (last_sgbc_helper_size_ != SGBCndofs_) {
            sgbc_helper_fields_ = initSGBCHelperFields(SGBCndofs_);
            last_sgbc_helper_size_ = SGBCndofs_;
        }

        mfem::Vector sgbcVec(6 * SGBCndofs_);
        sgbcVec = 0.0;

        for (auto& [tag, states] : sgbc_states_) {
            auto it = sgbc_wrapper_map_.find(tag);
            if (it == sgbc_wrapper_map_.end()) continue;

            auto w = it->second;
            for (const auto& state : states) {
                w->getSGBCFields(sgbc_sub_to_parent_ids_, state, sgbc_helper_fields_);

                const int global_L = state.global_pair.first;
                const int global_R = state.global_pair.second;
                const int sub_L = sgbc_sub_to_parent_ids_.Find(global_L);
                const int sub_R = sgbc_sub_to_parent_ids_.Find(global_R);

                if (sub_L == -1) continue;

                const int E_base_L = E*3*SGBCndofs_ + sub_L;
                const int H_base_L = H*3*SGBCndofs_ + sub_L;
                for (auto d : {X,Y,Z}) {
                    sgbcVec[E_base_L + d*SGBCndofs_] = sgbc_helper_fields_[E][d][sub_L];
                    sgbcVec[H_base_L + d*SGBCndofs_] = sgbc_helper_fields_[H][d][sub_L];
                }
                if (sub_R != -1){
                    const int E_base_R = E*3*SGBCndofs_ + sub_R;
                    const int H_base_R = H*3*SGBCndofs_ + sub_R;
                    for (auto d : {X,Y,Z}) {
                        sgbcVec[E_base_R + d*SGBCndofs_] = sgbc_helper_fields_[E][d][sub_R];
                        sgbcVec[H_base_R + d*SGBCndofs_] = sgbc_helper_fields_[H][d][sub_R];
                    }
                }
            }
        }

        mfem::Vector tempSGBC(sgbcVec.Size());
        tempSGBC.UseDevice(true);
        SGBCOperator_->Mult(sgbcVec, tempSGBC);

        for (int v = 0; v < sgbc_sub_to_parent_ids_.Size(); ++v) {
            int parent_idx = sgbc_sub_to_parent_ids_[v];
            const int E_out_base = E*3*ndofs + parent_idx;
            const int H_out_base = H*3*ndofs + parent_idx;
            const int E_temp_base = E*3*SGBCndofs_ + v;
            const int H_temp_base = H*3*SGBCndofs_ + v;
            for (int d : {X,Y,Z}) {
                out[E_out_base + d*ndofs] -= tempSGBC[E_temp_base + d*SGBCndofs_];
                out[H_out_base + d*ndofs] -= tempSGBC[H_temp_base + d*SGBCndofs_];
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
        }
        if (nbrDofs) {
            for (int d = X; d <= Z; ++d) {
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
        applyTFSFSourceToVector(GetTime(), ndofs, (srcmngr_.getGlobalTFSFSpace() ? srcmngr_.getGlobalTFSFSpace()->GetNDofs() : 0), s, true);
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