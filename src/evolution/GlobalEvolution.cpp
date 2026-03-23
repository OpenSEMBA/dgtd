#include "GlobalEvolution.h"
#include "solver/SolverExtension.h"
#include "components/SubMesher.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <unordered_set>
#ifdef SEMBA_DGTD_ENABLE_OPENMP
#include <omp.h>
#endif

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

    struct RawSGBCPair { NodePair pair; std::array<double,9> rot; };
    std::map<GeomTag, std::vector<RawSGBCPair>> raw_pairs;
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
                        std::array<double,9> id_rot = {1,0,0, 0,1,0, 0,0,1};
                        raw_pairs[bdr_attr].push_back({
                            NodePair(faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder(),
                                     faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder() + 1),
                            id_rot
                        });
                    } else {
                        auto twoElemSubMesh{ assembleInteriorFaceSubMesh(*mesh, *faceTrans, attMap) };
                        mfem::FiniteElementSpace subFES(&twoElemSubMesh, fec);
                        auto node_pair_local = buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES);

                        // Ensure consistent orientation: pair.first = inward side,
                        // pair.second = outward side (face normal direction).
                        // If face_ori < 0, Elem1 is on the outward side, so swap.
                        bool swap = buildFaceOrientation(*mesh, b) < 0.0;

                        // Compute face normal and build rotation matrix
                        auto face_normal = getNormal(*const_cast<mfem::FaceElementTransformations*>(faceTrans));
                        if (swap) face_normal *= -1.0; // ensure normal points from pair.first toward pair.second
                        std::array<double,9> face_rot{};
                        buildFaceRotationMatrix(face_normal, mesh->Dimension(), face_rot.data());

                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            if (swap) {
                                raw_pairs[bdr_attr].push_back({
                                    NodePair(node_pair_local.second[p], node_pair_local.first[p]),
                                    face_rot
                                });
                            } else {
                                raw_pairs[bdr_attr].push_back({
                                    NodePair(node_pair_local.first[p], node_pair_local.second[p]),
                                    face_rot
                                });
                            }
                        }
                    }
                }

                if (is_boundary) {
                    if (mesh->Dimension() == 1){
                        std::array<double,9> id_rot = {1,0,0, 0,1,0, 0,0,1};
                        if(faceTrans->Elem1No == 0) {
                            raw_pairs[bdr_attr].push_back({NodePair(0, -1), id_rot});
                        }
                        else {
                            raw_pairs[bdr_attr].push_back({
                                NodePair(faceTrans->Elem1No * (fec->GetOrder() + 1) + fec->GetOrder(), -1),
                                id_rot
                            });
                        }
                    } else {
                        auto elementSubMesh{ assembleBoundaryFaceSubMesh(*mesh, *faceTrans, attMap)};
                        mfem::FiniteElementSpace subFES(&elementSubMesh, fec);
                        auto node_pair_local = buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES);

                        auto face_normal = getNormal(*const_cast<mfem::FaceElementTransformations*>(faceTrans));
                        std::array<double,9> face_rot{};
                        buildFaceRotationMatrix(face_normal, mesh->Dimension(), face_rot.data());

                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            raw_pairs[bdr_attr].push_back({
                                NodePair(node_pair_local.first[p], -1),
                                face_rot
                            });
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

    // Diagnostic: print SGBC DOF pair info per rank
    for (const auto& [tag, pairs] : raw_pairs) {
        int n_interior = 0, n_boundary = 0;
        for (const auto& rp : pairs) {
            if (rp.pair.second != -1) ++n_interior; else ++n_boundary;
        }
        std::cout << "[SGBC Init] Rank " << Mpi::WorldRank()
                  << " | Tag " << tag
                  << " | " << pairs.size() << " DOF pairs"
                  << " (" << n_interior << " interior, " << n_boundary << " boundary)"
                  << " | Normal of first pair: ("
                  << pairs[0].rot[0] << ", " << pairs[0].rot[1] << ", " << pairs[0].rot[2] << ")"
                  << std::endl;
    }

    for (const auto& [tag, pairs] : raw_pairs){
        const auto& sbcps = model_.getSGBCProperties();
        for (auto p = 0; p < sbcps.size(); p++){
            for (auto t = 0; t < sbcps[p].geom_tags.size(); t++){
                if (tag == sbcps[p].geom_tags[t]){

                    auto wrap = SGBCWrapper::buildSGBCWrapper(sbcps[p]);
                    SGBCWrapper* wrap_ptr = wrap.get();

                    int state_size = wrap->getStateSize();
                    for(const auto& rp : pairs) {
                        SGBCState s;
                        s.global_pair = rp.pair;

                        s.element_pair.first = temp_dof_to_elem[rp.pair.first];
                        if (rp.pair.second != -1){
                            s.element_pair.second = temp_dof_to_elem[rp.pair.second];
                        }
                        else {
                            s.element_pair.second = -1;
                        }

                        if (s.element_pair.first == -1) {
                            std::cerr << "Error: Could not find Element ID for SGBC Node Pair: "
                                      << rp.pair.first << ", " << rp.pair.second << std::endl;
                        }

                        s.init(state_size);
                        std::memcpy(s.rot, rp.rot.data(), 9 * sizeof(double));
                        sgbc_states_[tag].push_back(s);
                    }

                    sgbcWrappers_.push_back(std::move(wrap));
                    sgbc_wrapper_map_[tag] = wrap_ptr;
                }
            }
        }
    }
    
    MPI_Barrier(fes_.GetParMesh()->GetComm());

    // Build flattened SGBC task list and per-thread wrapper clones for OpenMP.
    {
        for (const auto& [tag, states] : sgbc_states_) {
            for (size_t i = 0; i < states.size(); ++i) {
                sgbc_tasks_.push_back({tag, i});
            }
        }

#ifdef SEMBA_DGTD_ENABLE_OPENMP
        int max_threads = 1;
        #pragma omp parallel
        {
            #pragma omp single
            max_threads = omp_get_num_threads();
        }
        sgbc_omp_threads_ = std::min(max_threads, static_cast<int>(sgbc_tasks_.size()));
        if (sgbc_omp_threads_ > 1) {
            for (auto& [tag, _] : sgbc_wrapper_map_) {
                auto& pool = sgbc_thread_pool_[tag];
                pool.resize(sgbc_omp_threads_);
                // Thread 0 will use the original wrapper (no clone needed)
                for (int t = 1; t < sgbc_omp_threads_; ++t) {
                    pool[t] = sgbc_wrapper_map_[tag]->clone();
                }
            }
            if (Mpi::WorldRank() == 0) {
                std::cout << "[SGBC] OpenMP: using " << sgbc_omp_threads_
                          << " of " << max_threads << " available threads for "
                          << sgbc_tasks_.size() << " tasks" << std::endl;
            }
        } else if (Mpi::WorldRank() == 0 && !sgbc_tasks_.empty()) {
            std::cout << "[SGBC] OpenMP: 1 task, running serial" << std::endl;
        }
#endif
    }

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
        srcmngr_.initDirectPlanewaveEval();
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

GlobalEvolution::~GlobalEvolution() = default;

void GlobalEvolution::commitSGBCCheckpoint(double base_time, double dt)
{
    sgbc_step_base_time_ = base_time;
    sgbc_step_dt_ = dt;

    // First call: structurally clone the checkpoint map (allocate once)
    if (sgbc_states_checkpoint_.empty()) {
        for (const auto& [tag, states] : sgbc_states_) {
            auto& cp = sgbc_states_checkpoint_[tag];
            cp.resize(states.size());
            for (size_t i = 0; i < states.size(); ++i) {
                cp[i].global_pair  = states[i].global_pair;
                cp[i].element_pair = states[i].element_pair;
                cp[i].fields_state.SetSize(states[i].fields_state.Size());
            }
        }
    }
    // Fast path: only copy field data (no alloc/dealloc)
    for (const auto& [tag, states] : sgbc_states_) {
        auto& cp = sgbc_states_checkpoint_[tag];
        for (size_t i = 0; i < states.size(); ++i) {
            std::memcpy(cp[i].fields_state.GetData(),
                        states[i].fields_state.GetData(),
                        states[i].fields_state.Size() * sizeof(double));
        }
    }
}

void GlobalEvolution::finalizeSGBCStep(
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields)
{
    // Restore from checkpoint (memcpy only, sizes are pre-matched)
    for (auto& [tag, states] : sgbc_states_) {
        const auto& cp = sgbc_states_checkpoint_.at(tag);
        for (size_t i = 0; i < states.size(); ++i) {
            std::memcpy(states[i].fields_state.GetData(),
                        cp[i].fields_state.GetData(),
                        cp[i].fields_state.Size() * sizeof(double));
        }
    }

    // Advance SGBC for full dt with committed post-step fields
    double dt = sgbc_step_dt_;
#ifdef SEMBA_DGTD_ENABLE_OPENMP
    if (!sgbc_thread_pool_.empty()) {
        #pragma omp parallel for schedule(dynamic) num_threads(sgbc_omp_threads_)
        for (size_t ti = 0; ti < sgbc_tasks_.size(); ++ti) {
            const auto& task = sgbc_tasks_[ti];
            auto& state = sgbc_states_[task.tag][task.state_index];

            int tid = omp_get_thread_num();
            SGBCWrapper* w;
            if (tid == 0) {
                w = sgbc_wrapper_map_[task.tag];
            } else {
                w = sgbc_thread_pool_[task.tag][tid].get();
            }

            double sub_dt = w->getRecommendedDt();
            int nsteps = (sub_dt > 0.0 && sub_dt < dt)
                       ? static_cast<int>(std::ceil(dt / sub_dt))
                       : 1;
            double actual_sub_dt = dt / nsteps;

            w->loadState(state);
            for (int step = 0; step < nsteps; ++step) {
                w->updateFieldsWithGlobal(fields, state);
                w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
            }
            w->saveState(state);
        }
        for (auto& [tag, _] : sgbc_wrapper_map_) {
            _->setOldTime(sgbc_step_base_time_ + dt);
        }
    } else
#endif
    {
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
                for (int step = 0; step < nsteps; ++step) {
                    w->updateFieldsWithGlobal(fields, state);
                    w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
                }
                w->saveState(state);
            }
            w->setOldTime(sgbc_step_base_time_ + dt);
        }
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

    // Ensure global-layout planewave vector is allocated (first call only)
    if (tfsf_assembledFunc_.Size() != 6 * blockSize) {
        tfsf_assembledFunc_.SetSize(6 * blockSize);
        tfsf_assembledFunc_.UseDevice(true);
        tfsf_assembledFunc_ = 0.0;
    }
    // Vector is already zero: either freshly initialized above,
    // or cleaned after the previous AddMult call below.

    // Evaluate planewave on the TFSF submesh (with TF/SF masking and 0.5 scaling).
    // Fast path avoids ProjectCoefficient + deep copy of 6 GridFunctions.
    FieldGridFuncs func_storage;  // only used for fallback path
    const FieldGridFuncs* func_ptr;
    if (srcmngr_.hasDirectEval()) {
        srcmngr_.evalTimeVarFieldDirect(t_stage);
        func_ptr = &srcmngr_.getCachedTFSFFields();
    } else {
        func_storage = evalTimeVarFunction(t_stage, srcmngr_);
        func_ptr = &func_storage;
    }
    const auto& func = *func_ptr;

    // Scatter submesh planewave values into global-layout vector
    // Iterate DOFs in the outer loop to load each parent index only once
    for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v) {
        const int parent_dof = tfsf_sub_to_parent_ids_[v];
        for (int f : {E, H}) {
            for (int d : {X, Y, Z}) {
                tfsf_assembledFunc_[(f * 3 + d) * blockSize + parent_dof] = func[f][d][v];
            }
        }
    }

    // Apply TFSF operator and subtract from result: out -= TFSFOp * planewave
    TFSFOperator_->AddMult(tfsf_assembledFunc_, result_vector, -1.0);

    // Restore zeroed state for next call (only touch submesh entries)
    for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v) {
        const int parent_dof = tfsf_sub_to_parent_ids_[v];
        for (int comp = 0; comp < 6; ++comp) {
            tfsf_assembledFunc_[comp * blockSize + parent_dof] = 0.0;
        }
    }
}

void GlobalEvolution::Mult(const mfem::Vector& in, mfem::Vector& out) const
{
#ifdef SHOW_TIMER_INFORMATION
    mfem::StopWatch timerTotal, timerExchange, timerApplyA, timerTFSF, timerSGBC;
    timerTotal.Start();
#endif

    const auto ndofs     = fes_.GetNDofs();
    const auto nbrDofs   = fes_.num_face_nbr_dofs;
    const auto blockSize = ndofs + nbrDofs;

    // 1) SGBC sub-solve FIRST — independent of neighbor data, overlaps with
    //    other ranks' local work so they don't idle-wait during exchange.
#ifdef SHOW_TIMER_INFORMATION
    timerSGBC.Start();
#endif
    if (SGBCOperator_ && !sgbcWrappers_.empty()){
        if (sgbcVec_.Size() != 6 * blockSize) {
            sgbcVec_.SetSize(6 * blockSize);
            sgbcVec_.UseDevice(true);
            sgbcVec_ = 0.0;
        }

        double t_stage = GetTime();
        double dt_to_stage = t_stage - sgbc_step_base_time_;

        // Restore SGBC states from start-of-step checkpoint (memcpy only)
        for (auto& [tag, states] : sgbc_states_) {
            const auto& cp = sgbc_states_checkpoint_.at(tag);
            for (size_t i = 0; i < states.size(); ++i) {
                std::memcpy(states[i].fields_state.GetData(),
                            cp[i].fields_state.GetData(),
                            cp[i].fields_state.Size() * sizeof(double));
            }
        }

        // Advance SGBC to stage time using intermediate global fields as ghost data
        if (dt_to_stage > 1e-12) {
#ifdef SEMBA_DGTD_ENABLE_OPENMP
            if (!sgbc_thread_pool_.empty()) {
                #pragma omp parallel for schedule(dynamic) num_threads(sgbc_omp_threads_)
                for (size_t ti = 0; ti < sgbc_tasks_.size(); ++ti) {
                    const auto& task = sgbc_tasks_[ti];
                    auto& state = sgbc_states_[task.tag][task.state_index];

                    int tid = omp_get_thread_num();
                    SGBCWrapper* w;
                    if (tid == 0) {
                        w = sgbc_wrapper_map_.at(task.tag);
                    } else {
                        w = sgbc_thread_pool_.at(task.tag)[tid].get();
                    }

                    double sub_dt = w->getRecommendedDt();
                    int nsteps = (sub_dt > 0.0 && sub_dt < dt_to_stage)
                               ? static_cast<int>(std::ceil(dt_to_stage / sub_dt))
                               : 1;
                    double actual_sub_dt = dt_to_stage / nsteps;

                    w->loadState(state);
                    for (int step = 0; step < nsteps; ++step) {
                        w->updateFieldsWithGlobalVector(in, ndofs, state);
                        w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
                    }
                    w->saveState(state);
                }
            } else
#endif
            {
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
                        for (int step = 0; step < nsteps; ++step) {
                            w->updateFieldsWithGlobalVector(in, ndofs, state);
                            w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
                        }
                        w->saveState(state);
                    }
                }
            }
        }
    }
#ifdef SHOW_TIMER_INFORMATION
    timerSGBC.Stop();
#endif

    // 2) Exchange face neighbor data (MPI)
#ifdef SHOW_TIMER_INFORMATION
    timerExchange.Start();
#endif
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
#ifdef SHOW_TIMER_INFORMATION
    timerExchange.Stop();
#endif

    // 3) Assemble inNew_ (local + neighbor data)
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

    // 4) Apply globalOperator_ (DG flux) 
#ifdef SHOW_TIMER_INFORMATION
    timerApplyA.Start();
#endif
    if (out.Size() != 6 * ndofs) {
        out.SetSize(6 * ndofs);
        out.UseDevice(true);
    }
    globalOperator_->Mult(inNew_, out);
#ifdef SHOW_TIMER_INFORMATION
    timerApplyA.Stop();
#endif

    // 5) TFSF source injection
#ifdef SHOW_TIMER_INFORMATION
    timerTFSF.Start();
#endif
    applyTFSFSourceToVector(GetTime(), ndofs, nbrDofs, out);
#ifdef SHOW_TIMER_INFORMATION
    timerTFSF.Stop();
#endif

    // 6) SGBC flux injection (two-pass face coupling)
    //    Sub-solve was done in step 1; here we only inject interface values as flux.
    if (SGBCOperator_ && !sgbcWrappers_.empty()){
        if (sgbcVec_.Size() != 6 * blockSize) {
            sgbcVec_.SetSize(6 * blockSize);
            sgbcVec_.UseDevice(true);
            sgbcVec_ = 0.0;
        }
        for (auto& [tag, states] : sgbc_states_) {
            auto it = sgbc_wrapper_map_.find(tag);
            if (it == sgbc_wrapper_map_.end()) continue;
            auto w = it->second;
            int local_size = w->getLocalFieldSize();
            int idx_left = w->getLeftInterfaceIndex();
            const auto left_bdr =
                w->getProperties().sgbc_bdr_info.first.bdrCond;

            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                const double* R = state.rot;
                // Elem1 slot: global field value (u_L) — already in global frame
                for (int d = X; d <= Z; ++d) {
                    sgbcVec_[d * blockSize + gl] = in[d * ndofs + gl];
                    sgbcVec_[(3 + d) * blockSize + gl] = in[(3 + d) * ndofs + gl];
                }
                // Elem2 slot: what Elem1 sees from the slab's left side
                if (gr != -1) {
                    if (left_bdr == BdrCond::PEC) {
                        // PEC on the left: Elem1 sees a PEC wall.
                        for (int d = X; d <= Z; ++d) {
                            sgbcVec_[d * blockSize + gr]       = -in[d * ndofs + gl];
                            sgbcVec_[(3 + d) * blockSize + gr] =  in[(3 + d) * ndofs + gl];
                        }
                    } else if (left_bdr == BdrCond::PMC) {
                        // PMC on the left: Elem1 sees a PMC wall.
                        for (int d = X; d <= Z; ++d) {
                            sgbcVec_[d * blockSize + gr]       =  in[d * ndofs + gl];
                            sgbcVec_[(3 + d) * blockSize + gr] = -in[(3 + d) * ndofs + gl];
                        }
                    } else if (left_bdr == BdrCond::SMA) {
                        // SMA on the left: Elem1 sees a matched load.
                        for (int d = X; d <= Z; ++d) {
                            sgbcVec_[d * blockSize + gr]       = in[d * ndofs + gl];
                            sgbcVec_[(3 + d) * blockSize + gr] = in[(3 + d) * ndofs + gl];
                        }
                    } else {
                        // Transparent: use slab left-interface value — rotate local→global
                        double el[3], hl[3];
                        for (int d = 0; d < 3; ++d) {
                            el[d] = state.fields_state[d * local_size + idx_left];
                            hl[d] = state.fields_state[(3 + d) * local_size + idx_left];
                        }
                        for (int d = 0; d < 3; ++d) {
                            sgbcVec_[d * blockSize + gr]       = R[d]*el[0] + R[3+d]*el[1] + R[6+d]*el[2];
                            sgbcVec_[(3 + d) * blockSize + gr] = R[d]*hl[0] + R[3+d]*hl[1] + R[6+d]*hl[2];
                        }
                    }
                }
            }
        }
        // Targeted SpMV: only compute rows for Elem1 DOFs, add directly to out.
        {
            const double* sv = sgbcVec_.HostRead();
            double* op = out.HostReadWrite();
            for (int comp = 0; comp < 6; ++comp) {
                for (int dof : sgbc_elem1_dofs_) {
                    const int row = comp * ndofs + dof;
                    const int nnz = SGBCOperator_->RowSize(row);
                    const int* cols = SGBCOperator_->GetRowColumns(row);
                    const double* vals = SGBCOperator_->GetRowEntries(row);
                    double sum = 0.0;
                    for (int j = 0; j < nnz; ++j) {
                        sum += vals[j] * sv[cols[j]];
                    }
                    op[row] += sum;
                }
            }
        }

        // Clear only the entries that were written in Pass 1
        for (auto& [tag, states] : sgbc_states_) {
            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                for (int d = X; d <= Z; ++d) {
                    sgbcVec_[d * blockSize + gl] = 0.0;
                    sgbcVec_[(3 + d) * blockSize + gl] = 0.0;
                    if (gr != -1) {
                        sgbcVec_[d * blockSize + gr] = 0.0;
                        sgbcVec_[(3 + d) * blockSize + gr] = 0.0;
                    }
                }
            }
        }

        // --- Pass 2: Elem2 (R) flux ---
        // Input: u_slab_right at gl (Elem1), u_R at gr (Elem2)
        // Keep only Elem2 DOF contributions from the result
        for (auto& [tag, states] : sgbc_states_) {
            auto it = sgbc_wrapper_map_.find(tag);
            if (it == sgbc_wrapper_map_.end()) continue;
            auto w = it->second;
            int local_size = w->getLocalFieldSize();
            int idx_right = w->getRightInterfaceIndex();
            const auto right_bdr =
                w->getProperties().sgbc_bdr_info.second.bdrCond;

            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                if (gr == -1) continue;
                const double* R = state.rot;

                if (right_bdr == BdrCond::PEC) {
                    // PEC on the right: Elem2 sees a PEC wall.
                    // Mirror state: E_ext = -E_int, H_ext = +H_int
                    for (int d = X; d <= Z; ++d) {
                        sgbcVec_[d * blockSize + gl]       = -in[d * ndofs + gr];
                        sgbcVec_[(3 + d) * blockSize + gl] =  in[(3 + d) * ndofs + gr];
                    }
                } else if (right_bdr == BdrCond::PMC) {
                    // PMC on the right: Elem2 sees a PMC wall.
                    // Mirror state: E_ext = +E_int, H_ext = -H_int
                    for (int d = X; d <= Z; ++d) {
                        sgbcVec_[d * blockSize + gl]       =  in[d * ndofs + gr];
                        sgbcVec_[(3 + d) * blockSize + gl] = -in[(3 + d) * ndofs + gr];
                    }
                } else if (right_bdr == BdrCond::SMA) {
                    // SMA (absorbing) on the right: Elem2 sees a matched load.
                    // Exterior state = interior state → zero jump → no reflection.
                    for (int d = X; d <= Z; ++d) {
                        sgbcVec_[d * blockSize + gl]       = in[d * ndofs + gr];
                        sgbcVec_[(3 + d) * blockSize + gl] = in[(3 + d) * ndofs + gr];
                    }
                } else {
                    // Transparent: use slab right-interface value — rotate local→global
                    double el[3], hl[3];
                    for (int d = 0; d < 3; ++d) {
                        el[d] = state.fields_state[d * local_size + idx_right];
                        hl[d] = state.fields_state[(3 + d) * local_size + idx_right];
                    }
                    for (int d = 0; d < 3; ++d) {
                        sgbcVec_[d * blockSize + gl]       = R[d]*el[0] + R[3+d]*el[1] + R[6+d]*el[2];
                        sgbcVec_[(3 + d) * blockSize + gl] = R[d]*hl[0] + R[3+d]*hl[1] + R[6+d]*hl[2];
                    }
                }
                // Elem2 slot: global field value (u_R) — already in global frame
                for (int d = X; d <= Z; ++d) {
                    sgbcVec_[d * blockSize + gr] = in[d * ndofs + gr];
                    sgbcVec_[(3 + d) * blockSize + gr] = in[(3 + d) * ndofs + gr];
                }
            }
        }
        // Targeted SpMV: only compute rows for Elem2 DOFs, add directly to out.
        {
            const double* sv = sgbcVec_.HostRead();
            double* op = out.HostReadWrite();
            for (int comp = 0; comp < 6; ++comp) {
                for (int dof : sgbc_elem2_dofs_) {
                    const int row = comp * ndofs + dof;
                    const int nnz = SGBCOperator_->RowSize(row);
                    const int* cols = SGBCOperator_->GetRowColumns(row);
                    const double* vals = SGBCOperator_->GetRowEntries(row);
                    double sum = 0.0;
                    for (int j = 0; j < nnz; ++j) {
                        sum += vals[j] * sv[cols[j]];
                    }
                    op[row] += sum;
                }
            }
        }

        // Clear entries written in Pass 2 so sgbcVec_ stays zeroed for next call
        for (auto& [tag, states] : sgbc_states_) {
            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                if (gr == -1) continue;
                for (int d = X; d <= Z; ++d) {
                    sgbcVec_[d * blockSize + gl] = 0.0;
                    sgbcVec_[(3 + d) * blockSize + gl] = 0.0;
                    sgbcVec_[d * blockSize + gr] = 0.0;
                    sgbcVec_[(3 + d) * blockSize + gr] = 0.0;
                }
            }
        }
    }

#ifdef SHOW_TIMER_INFORMATION
    timerTotal.Stop();

    // Skip timer output for SGBC sub-solver Mult() calls.
    // In multi-rank runs the MPI collectives would deadlock because only
    // SGBC ranks enter the sub-solver path.
    if (!opts_.is_sgbc_solver) {
        static double acc_total = 0, acc_sgbc = 0, acc_exchange = 0, acc_applyA = 0, acc_tfsf = 0;
        static int acc_count = 0;
        static double lastPrintTime = -1.0;
        constexpr double printInterval = 0.05;
        constexpr int printCadence = 50;

        acc_total    += timerTotal.RealTime()    * 1000.0;
        acc_sgbc     += timerSGBC.RealTime()     * 1000.0;
        acc_exchange += timerExchange.RealTime() * 1000.0;
        acc_applyA   += timerApplyA.RealTime()   * 1000.0;
        acc_tfsf     += timerTFSF.RealTime()     * 1000.0;
        acc_count++;

        const double currentTime = GetTime();
        bool want_print = (acc_count >= printCadence &&
            (lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval));

        const int P = Mpi::WorldSize();

        if (P > 1) {
            // Multi-rank: synchronize print decision across all ranks.
            int local_want = want_print ? 1 : 0;
            int global_want = 0;
            MPI_Allreduce(&local_want, &global_want, 1, MPI_INT, MPI_MAX,
                          fes_.GetParMesh()->GetComm());
            want_print = (global_want != 0);
        }

        if (want_print)
        {
            double local_avg[5] = {
                acc_total    / acc_count,
                acc_sgbc     / acc_count,
                acc_exchange / acc_count,
                acc_applyA   / acc_count,
                acc_tfsf     / acc_count
            };

            std::vector<double> all_avg(5 * P);
            if (P > 1) {
                MPI_Gather(local_avg, 5, MPI_DOUBLE,
                           all_avg.data(), 5, MPI_DOUBLE,
                           0, fes_.GetParMesh()->GetComm());
            } else {
                std::copy(local_avg, local_avg + 5, all_avg.data());
            }

            if (Mpi::WorldRank() == 0) {
                std::cout << "\n[Mult timing] t=" << std::fixed << std::setprecision(4) << currentTime
                          << "  (avg of " << acc_count << " calls, ms/call)\n"
                          << std::defaultfloat;
                std::cout << "  Rank  | total    sgbc   exchange  operator   tfsf  | bottleneck\n"
                          << "  ------+----------------------------------------------------+-----------\n";
                for (int r = 0; r < P; ++r) {
                    double t_total    = all_avg[r*5 + 0];
                    double t_sgbc     = all_avg[r*5 + 1];
                    double t_exchange = all_avg[r*5 + 2];
                    double t_applyA   = all_avg[r*5 + 3];
                    double t_tfsf     = all_avg[r*5 + 4];
                    const char* bottleneck = "compute";
                    if (t_exchange > 0.7 * t_total) bottleneck = "MPI wait";
                    else if (t_sgbc > 0.5 * t_total) bottleneck = "SGBC";

                    char buf[256];
                    std::snprintf(buf, sizeof(buf),
                        "  %4d  | %7.2f  %6.2f  %8.2f  %8.2f  %6.3f  | %s",
                        r, t_total, t_sgbc, t_exchange, t_applyA, t_tfsf, bottleneck);
                    std::cout << buf << '\n';
                }
                std::cout << std::flush;
            }

            acc_total = acc_sgbc = acc_exchange = acc_applyA = acc_tfsf = 0;
            acc_count = 0;
            lastPrintTime = currentTime;
        }
    }
#endif
}

void GlobalEvolution::ImplicitSolve(const double dt,
                                    const mfem::Vector& x,
                                    mfem::Vector& k)
{
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

    // --- Compute RHS = A*x + source ---
    // For serial meshes (nbrDofs==0), globalOperator_ is square n×n and the
    // copy-through-GridFunction round-trip is identity. Use direct SpMV.
    if (nbrDofs == 0) {
        globalOperator_->Mult(x, implicit_rhs_);
    } else {
        // Parallel path: need to exchange face neighbor data
        for (int d = X; d <= Z; ++d) {
            std::memcpy(eOld_[d].GetData(), x.GetData() + d * ndofs, ndofs * sizeof(double));
            std::memcpy(hOld_[d].GetData(), x.GetData() + (3 + d) * ndofs, ndofs * sizeof(double));
        }
        for (int d = X; d <= Z; ++d) {
            eOld_[d].ExchangeFaceNbrData();
            hOld_[d].ExchangeFaceNbrData();
        }
        for (int d = X; d <= Z; ++d) {
            implicit_inNew_.SetVector(eOld_[d],       d      * blockSize);
            implicit_inNew_.SetVector(hOld_[d],  (3 + d)     * blockSize);
        }
        for (int d = X; d <= Z; ++d) {
            mfem::Vector &eNbr = eOld_[d].FaceNbrData();
            mfem::Vector &hNbr = hOld_[d].FaceNbrData();
            implicit_inNew_.SetVector(eNbr,      d      * blockSize + ndofs);
            implicit_inNew_.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
        }
        globalOperator_->Mult(implicit_inNew_, implicit_rhs_);
    }

    // Add TFSF source contribution (skip if no active TFSF sources)
    if (!tfsf_cutoff_reached_ && TFSFOperator_ && !tfsfSourceIndices_.empty()) {
        implicit_src_ = 0.0;
        applyTFSFSourceToVector(GetTime(), ndofs, nbrDofs, implicit_src_);
        implicit_rhs_ += implicit_src_;
    }

    // --- Solve (I - dt*A) k = rhs ---
    // For small serial systems (SGBC sub-solver), use cached dense LU
    // factorization. One-time cost: form A as dense + LU factor.
    // Per-solve cost: O(n^2) back-substitution instead of ~30 GMRES iterations.
    if (nbrDofs == 0 && n <= dense_solve_threshold_) {
        // One-time: form dense copy of globalOperator_
        if (!dense_A_formed_) {
            globalOperator_->ToDenseMatrix(A_dense_);
            J_dense_.SetSize(n);
            dense_A_formed_ = true;
            cached_dt_ = -1.0;  // force re-factorization
        }

        // Re-factorize when dt changes (different sub-step sizes)
        if (dt != cached_dt_) {
            // J = I - dt*A
            J_dense_ = A_dense_;
            J_dense_ *= -dt;
            for (int i = 0; i < n; i++) {
                J_dense_(i, i) += 1.0;
            }
            J_inv_ = std::make_unique<mfem::DenseMatrixInverse>(J_dense_);
            cached_dt_ = dt;
        }

        if (k.Size() != n) {
            k.SetSize(n);
        }
        J_inv_->Mult(implicit_rhs_, k);
    } else {
        // Large/parallel system: fall back to GMRES
        auto applyA_parallel = [&](const mfem::Vector& u, mfem::Vector& Au)
        {
            for (int d = X; d <= Z; ++d) {
                std::memcpy(eOld_[d].GetData(), u.GetData() + d * ndofs, ndofs * sizeof(double));
                std::memcpy(hOld_[d].GetData(), u.GetData() + (3 + d) * ndofs, ndofs * sizeof(double));
            }
            for (int d = X; d <= Z; ++d) {
                eOld_[d].ExchangeFaceNbrData();
                hOld_[d].ExchangeFaceNbrData();
            }
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
            globalOperator_->Mult(implicit_inNew_, Au);
        };

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
        } J(n, std::function<void(const mfem::Vector&, mfem::Vector&)>(applyA_parallel), dt);

        mfem::GMRESSolver gmres;
        gmres.SetOperator(J);
        gmres.SetRelTol(1e-8);
        gmres.SetMaxIter(400);
        gmres.SetKDim(60);
        gmres.SetPrintLevel(0);

        if (k.Size() != n) {
            k.SetSize(n); k.UseDevice(true); k = implicit_rhs_;
        }
        gmres.Mult(implicit_rhs_, k);
    }
}

}