#include "GlobalEvolution.h"
#include "solver/SolverExtension.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <unordered_set>

namespace maxwell {

SGBCHelperFields initSGBCHelperFields(const int size)
{
    SGBCHelperFields res;
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

    MPI_Barrier(fes_.GetParMesh()->GetComm());

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
    
    MPI_Barrier(fes_.GetParMesh()->GetComm());

    // Pre-compute unique element DOF sets for SGBC flux gathering.
    // Instead of copying only face DOFs from sgbcOut_ to out, we need ALL DOFs
    // of elements touching SGBC faces, because the DG face integral distributes
    // flux to all test functions in the element (not just the face node).
    {
        std::unordered_set<int> elem1_set, elem2_set;
        for (const auto& [tag, states] : sgbc_states_) {
            for (const auto& s : states) {
                elem1_set.insert(s.element_pair.first);
                if (s.element_pair.second != -1) {
                    elem2_set.insert(s.element_pair.second);
                }
            }
        }
        const mfem::Table& elem_dof_table = fes_.GetElementToDofTable();
        auto collectDofs = [&](const std::unordered_set<int>& elems, std::vector<int>& out_dofs) {
            std::unordered_set<int> dof_set;
            for (int el : elems) {
                const int* dofs = elem_dof_table.GetRow(el);
                const int n = elem_dof_table.RowSize(el);
                for (int i = 0; i < n; ++i) {
                    dof_set.insert(dofs[i]);
                }
            }
            out_dofs.assign(dof_set.begin(), dof_set.end());
            std::sort(out_dofs.begin(), out_dofs.end());
        };
        collectDofs(elem1_set, sgbc_elem1_dofs_);
        collectDofs(elem2_set, sgbc_elem2_dofs_);
    }

    globalOperator_ = std::make_unique<mfem::SparseMatrix>(
        numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs(), 
        numberOfFieldComponents * numberOfMaxDimensions * (fes_.GetNDofs() + fes_.num_face_nbr_dofs)
    );

    Probes probes; 

    // Keep TFSF submesh infrastructure for planewave source evaluation
    if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {
        srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
        auto src_sm = static_cast<mfem::SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
        mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), tfsf_sub_to_parent_ids_);
    }

    // Build all operators on the global mesh
    ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
    DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, fes_);
    globalOperator_ = dgops.buildGlobalOperator();

    if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {
        TFSFOperator_ = dgops.buildTFSFGlobalOperator();
    }
    if (model_.getSGBCToMarker().find(BdrCond::SGBC) != model_.getSGBCToMarker().end()) {
        SGBCOperator_ = dgops.buildSGBCGlobalOperator();
    }

    // --- Performance: cache which sources are TotalField ---
    {
        int idx = 0;
        for (auto& source : srcmngr_.sources) {
            if (source && dynamic_cast<TotalField*>(source.get())) {
                tfsfSourceIndices_.push_back(idx);
            }
            ++idx;
        }
    }
}

const mfem::Vector buildSingleVectorTFSFFunc(const FieldGridFuncs& func)
{
    const int size = func[E][X].Size();
    mfem::Vector res(6 * size);
    res.UseDevice(true);
    for (int f : { E, H }) {
        for (int d : { X, Y, Z }) {
            std::memcpy(res.GetData() + (f * 3 + d) * size,
                        func[f][d].GetData(),
                        size * sizeof(double));
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

void GlobalEvolution::commitSGBCCheckpoint(double base_time, double dt)
{
    sgbc_step_base_time_ = base_time;
    sgbc_step_dt_ = dt;

    // Deep copy all SGBC states (each SGBCState contains an mfem::Vector)
    sgbc_states_checkpoint_.clear();
    for (const auto& [tag, states] : sgbc_states_) {
        auto& cp = sgbc_states_checkpoint_[tag];
        cp.reserve(states.size());
        for (const auto& s : states) {
            SGBCState copy;
            copy.global_pair = s.global_pair;
            copy.element_pair = s.element_pair;
            copy.fields_state = s.fields_state;  // mfem::Vector deep copy
            cp.push_back(std::move(copy));
        }
    }
}

void GlobalEvolution::finalizeSGBCStep(
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields)
{
    // Restore from checkpoint
    for (auto& [tag, states] : sgbc_states_) {
        const auto& cp = sgbc_states_checkpoint_.at(tag);
        for (size_t i = 0; i < states.size(); ++i) {
            states[i].fields_state = cp[i].fields_state;
        }
    }

    // Advance SGBC for full dt with committed post-step fields
    double dt = sgbc_step_dt_;
    for (auto& [tag, states] : sgbc_states_) {
        auto it = sgbc_wrapper_map_.find(tag);
        if (it == sgbc_wrapper_map_.end()) continue;

        auto w = it->second;
        double sub_dt = w->getRecommendedDt();
        int nsteps = (sub_dt > 0.0 && sub_dt < dt)
                   ? static_cast<int>(std::ceil(dt / sub_dt))
                   : 1;
        double actual_sub_dt = dt / nsteps;

        for (auto& state : states) {
            w->loadState(state);
            w->updateFieldsWithGlobal(fields, state);
            for (int step = 0; step < nsteps; ++step) {
                w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
            }
            w->saveState(state);
        }
        w->setOldTime(sgbc_step_base_time_ + dt);
    }
}

void GlobalEvolution::applyTFSFSourceToVector(double t_stage, int ndofs, int nbrDofs,
                                              mfem::Vector& result_vector) const
{
    if (tfsf_cutoff_reached_ || !TFSFOperator_) return;
    if (tfsfSourceIndices_.empty()) return;

    auto *tfsf_space = srcmngr_.getGlobalTFSFSpace();
    if (!tfsf_space) return;

    const double t_off = opts_.tfsf_cutoff_time;
    if (t_off != 0.0 && t_stage >= t_off) {
        tfsf_cutoff_reached_ = true;
        return;
    }

    const int blockSize = ndofs + nbrDofs;

    // Ensure global-layout planewave vector is allocated
    if (tfsf_assembledFunc_.Size() != 6 * blockSize) {
        tfsf_assembledFunc_.SetSize(6 * blockSize);
        tfsf_assembledFunc_.UseDevice(true);
    }
    tfsf_assembledFunc_ = 0.0;

    // Evaluate planewave on the TFSF submesh (with TF/SF masking and 0.5 scaling)
    auto func = evalTimeVarFunction(t_stage, srcmngr_);

    // Scatter submesh planewave values into global-layout vector
    for (int f : {E, H}) {
        for (int d : {X, Y, Z}) {
            for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v) {
                tfsf_assembledFunc_[(f * 3 + d) * blockSize + tfsf_sub_to_parent_ids_[v]] = func[f][d][v];
            }
        }
    }

    // Apply TFSF operator and subtract from result: out -= TFSFOp * planewave
    TFSFOperator_->AddMult(tfsf_assembledFunc_, result_vector, -1.0);
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
        std::memcpy(eOld_[d].GetData(), in.GetData() + d * ndofs, ndofs * sizeof(double));
        std::memcpy(hOld_[d].GetData(), in.GetData() + (3 + d) * ndofs, ndofs * sizeof(double));
    }
#endif
    for (int d = X; d <= Z; ++d) {
        eOld_[d].ExchangeFaceNbrData();
        hOld_[d].ExchangeFaceNbrData();
    }
    timerExchange.Stop();

    timerAssembleInNew.Start();
    if (inNew_.Size() != 6 * blockSize) {
        inNew_.SetSize(6 * blockSize);
        inNew_.UseDevice(true);
    }

#ifdef SEMBA_DGTD_ENABLE_CUDA
    load_eh_to_innew_gpu(in, inNew_, ndofs, nbrDofs);
    load_nbr_to_innew_gpu(eOld_, hOld_, inNew_, ndofs, nbrDofs);
#else
    for (int d = X; d <= Z; ++d) {
        inNew_.SetVector(eOld_[d],       d      * blockSize);
        inNew_.SetVector(hOld_[d],  (3 + d)     * blockSize);
    }
    if (nbrDofs) {
        for (int d = X; d <= Z; ++d) {
            mfem::Vector &eNbr = eOld_[d].FaceNbrData();
            mfem::Vector &hNbr = hOld_[d].FaceNbrData();
            inNew_.SetVector(eNbr,      d      * blockSize + ndofs);
            inNew_.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
        }
    }
#endif
    timerAssembleInNew.Stop();

    timerApplyA.Start();
    if (out.Size() != 6 * ndofs) {
        out.SetSize(6 * ndofs);
        out.UseDevice(true);
    }
    globalOperator_->Mult(inNew_, out);
    timerApplyA.Stop();

    timerTFSF.Start();
    applyTFSFSourceToVector(GetTime(), ndofs, nbrDofs, out);
    timerTFSF.Stop();

    // 5) SGBC Coupling via Monolithic IMEX: two-pass face flux injection
    //
    // The global operator now EXCLUDES SGBC faces from its interior face flux.
    // We compute the correct two-interface flux here:
    //   Pass 1: Elem1 (L) sees itself + slab_left as neighbor -> keep Elem1 DOF contributions
    //   Pass 2: Elem2 (R) sees slab_right as neighbor + itself -> keep Elem2 DOF contributions
    timerSGBC.Start();

    if (SGBCOperator_ && !sgbcWrappers_.empty()){
        if (sgbcVec_.Size() != 6 * blockSize) {
            sgbcVec_.SetSize(6 * blockSize);
            sgbcVec_.UseDevice(true);
        }
        if (sgbcOut_.Size() != 6 * ndofs) {
            sgbcOut_.SetSize(6 * ndofs);
            sgbcOut_.UseDevice(true);
        }

        double t_stage = GetTime();
        double dt_to_stage = t_stage - sgbc_step_base_time_;

        // Restore SGBC states from start-of-step checkpoint
        for (auto& [tag, states] : sgbc_states_) {
            const auto& cp = sgbc_states_checkpoint_.at(tag);
            for (size_t i = 0; i < states.size(); ++i) {
                states[i].fields_state = cp[i].fields_state;
            }
        }

        // Advance SGBC to stage time using intermediate global fields as ghost data
        if (dt_to_stage > 1e-12) {
            for (auto& [tag, states] : sgbc_states_) {
                auto it = sgbc_wrapper_map_.find(tag);
                if (it == sgbc_wrapper_map_.end()) continue;

                auto w = it->second;
                double sub_dt = w->getRecommendedDt();
                int nsteps = (sub_dt > 0.0 && sub_dt < dt_to_stage)
                           ? static_cast<int>(std::ceil(dt_to_stage / sub_dt))
                           : 1;
                double actual_sub_dt = dt_to_stage / nsteps;

                for (auto& state : states) {
                    w->loadState(state);
                    w->updateFieldsWithGlobalVector(in, ndofs, state);
                    for (int step = 0; step < nsteps; ++step) {
                        w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
                    }
                    w->saveState(state);
                }
            }
        }

        // --- Pass 1: Elem1 (L) flux ---
        // Input: u_L at gl (Elem1), u_slab_left at gr (Elem2)
        // Keep only Elem1 DOF contributions from the result
        sgbcVec_ = 0.0;
        for (auto& [tag, states] : sgbc_states_) {
            auto it = sgbc_wrapper_map_.find(tag);
            if (it == sgbc_wrapper_map_.end()) continue;
            auto w = it->second;
            int local_size = w->getLocalFieldSize();
            int idx_left = w->getLeftInterfaceIndex();

            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                for (int d = X; d <= Z; ++d) {
                    // Elem1 slot: global field value (u_L)
                    sgbcVec_[d * blockSize + gl] = in[d * ndofs + gl];
                    sgbcVec_[(3 + d) * blockSize + gl] = in[(3 + d) * ndofs + gl];
                    // Elem2 slot: slab left-interface value
                    if (gr != -1) {
                        sgbcVec_[d * blockSize + gr] = state.fields_state[d * local_size + idx_left];
                        sgbcVec_[(3 + d) * blockSize + gr] = state.fields_state[(3 + d) * local_size + idx_left];
                    }
                }
            }
        }
        sgbcOut_ = 0.0;
        SGBCOperator_->Mult(sgbcVec_, sgbcOut_);
        // Add ALL Elem1 DOF contributions to out (not just face DOFs).
        // The DG face integral distributes flux to all test functions in
        // the element, so we must copy the full element contribution.
        for (int comp = 0; comp < 6; ++comp) {
            for (int dof : sgbc_elem1_dofs_) {
                out[comp * ndofs + dof] += sgbcOut_[comp * ndofs + dof];
            }
        }

        // --- Pass 2: Elem2 (R) flux ---
        // Input: u_slab_right at gl (Elem1), u_R at gr (Elem2)
        // Keep only Elem2 DOF contributions from the result
        sgbcVec_ = 0.0;
        for (auto& [tag, states] : sgbc_states_) {
            auto it = sgbc_wrapper_map_.find(tag);
            if (it == sgbc_wrapper_map_.end()) continue;
            auto w = it->second;
            int local_size = w->getLocalFieldSize();
            int idx_right = w->getRightInterfaceIndex();

            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                if (gr == -1) continue;
                for (int d = X; d <= Z; ++d) {
                    // Elem1 slot: slab right-interface value
                    sgbcVec_[d * blockSize + gl] = state.fields_state[d * local_size + idx_right];
                    sgbcVec_[(3 + d) * blockSize + gl] = state.fields_state[(3 + d) * local_size + idx_right];
                    // Elem2 slot: global field value (u_R)
                    sgbcVec_[d * blockSize + gr] = in[d * ndofs + gr];
                    sgbcVec_[(3 + d) * blockSize + gr] = in[(3 + d) * ndofs + gr];
                }
            }
        }
        sgbcOut_ = 0.0;
        SGBCOperator_->Mult(sgbcVec_, sgbcOut_);
        // Add ALL Elem2 DOF contributions to out (not just face DOFs).
        for (int comp = 0; comp < 6; ++comp) {
            for (int dof : sgbc_elem2_dofs_) {
                out[comp * ndofs + dof] += sgbcOut_[comp * ndofs + dof];
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

    // Reuse work vectors across calls (avoid repeated allocation)
    if (!implicit_work_initialized_ || implicit_inNew_.Size() != 6 * blockSize) {
        implicit_inNew_.SetSize(6 * blockSize);
        implicit_inNew_.UseDevice(true);
        implicit_rhs_.SetSize(n);
        implicit_rhs_.UseDevice(true);
        implicit_src_.SetSize(n);
        implicit_src_.UseDevice(true);
        implicit_work_initialized_ = true;
    }

    auto applyA_noPrint = [&](const mfem::Vector& u, mfem::Vector& Au)
    {
#ifdef SEMBA_DGTD_ENABLE_CUDA
        load_in_to_eh_gpu(u, eOld_, hOld_, ndofs);
#else
        for (int d = X; d <= Z; ++d) {
            std::memcpy(eOld_[d].GetData(), u.GetData() + d * ndofs, ndofs * sizeof(double));
            std::memcpy(hOld_[d].GetData(), u.GetData() + (3 + d) * ndofs, ndofs * sizeof(double));
        }
#endif
        for (int d = X; d <= Z; ++d) {
            eOld_[d].ExchangeFaceNbrData();
            hOld_[d].ExchangeFaceNbrData();
        }

#ifdef SEMBA_DGTD_ENABLE_CUDA
        load_eh_to_innew_gpu(u, implicit_inNew_, ndofs, nbrDofs);
        load_nbr_to_innew_gpu(eOld_, hOld_, implicit_inNew_, ndofs, nbrDofs);
#else
        for (int d = X; d <= Z; ++d) {
            implicit_inNew_.SetVector(eOld_[d],       d      * blockSize);
            implicit_inNew_.SetVector(hOld_[d],  (3 + d)     * blockSize);
        }
        if (nbrDofs) {
            for (int d = X; d <= Z; ++d) {
                mfem::Vector &eNbr = eOld_[d].FaceNbrData();
                mfem::Vector &hNbr = hOld_[d].FaceNbrData();
                implicit_inNew_.SetVector(eNbr,      d      * blockSize + ndofs);
                implicit_inNew_.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
            }
        }
#endif

        globalOperator_->Mult(implicit_inNew_, Au);
    };

    timerRHS.Start();
    timerAx.Start();
    applyA_noPrint(x, implicit_rhs_);
    timerAx.Stop();

    timerSrc.Start();
    {
        implicit_src_ = 0.0;
        applyTFSFSourceToVector(GetTime(), ndofs, nbrDofs, implicit_src_);
        implicit_rhs_ += implicit_src_;
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

    if (k.Size() != n) {
        k.SetSize(n); k.UseDevice(true); k = implicit_rhs_;
    }

    timerSolve.Start();
    gmres.Mult(implicit_rhs_, k);
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