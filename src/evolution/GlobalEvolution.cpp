#include "GlobalEvolution.h"
#include "components/PMLOperatorAudit.h"
#include "MaxwellEvolutionMethods.h"
#include "components/PMLDGHelpers.h"
#include "components/PMLSignTest.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include <fstream>
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

static std::vector<int> collectRCSSurfaceTags(const Probes& probes)
{
    std::unordered_set<int> unique_tags;
    for (const auto& p : probes.rcsSurfaceProbes) {
        for (const int t : p.tags) {
            unique_tags.insert(t);
        }
    }
    std::vector<int> tags(unique_tags.begin(), unique_tags.end());
    std::sort(tags.begin(), tags.end());
    return tags;
}

static void logTFSFPartitionFaceCoverage(mfem::ParMesh& mesh,
                                         const mfem::Array<int>& marker)
{
    if (marker.Size() == 0) return;
    bool any_tag = false;
    for (int i = 0; i < marker.Size(); ++i) {
        if (marker[i] != 0) { any_tag = true; break; }
    }
    if (!any_tag) return;

    int local_ibfi = 0;
    int local_bfi = 0;
    for (int b = 0; b < mesh.GetNBE(); ++b) {
        const int attr = mesh.GetBdrAttribute(b) - 1;
        if (attr < 0 || attr >= marker.Size() || marker[attr] == 0) continue;
        if (mesh.GetInternalBdrFaceTransformations(b) != nullptr) {
            ++local_ibfi;
        } else {
            ++local_bfi;
        }
    }

    const MPI_Comm comm = mesh.GetComm();
    int global_ibfi = 0;
    int global_bfi = 0;
    MPI_Allreduce(&local_ibfi, &global_ibfi, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&local_bfi, &global_bfi, 1, MPI_INT, MPI_SUM, comm);

    if (Mpi::WorldRank() == 0) {
        std::cout << "[TFSF] face coverage across ranks: IBFI=" << global_ibfi
                  << " BFI(partition-boundary)=" << global_bfi;
        if (global_bfi > 0) {
            std::cout << "  (BFI path active for MPI partition faces)";
        }
        std::cout << std::endl;
    }
}

GlobalEvolution::GlobalEvolution(
    mfem::ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options, const Probes& probes, double final_time) :
    mfem::TimeDependentOperator([&]() {
        const int ndofs = fes.GetNDofs();
        const int n_aux = model.hasPML() && model.getPMLAuxLayout()
            ? model.getPMLAuxLayout()->nAux()
            : 0;
        return numberOfFieldComponents * numberOfMaxDimensions * ndofs + n_aux;
    }()),
    total_state_size_(Height()),
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

                        // Compute per-DOF normals on curved faces.
                        // Find local DOF indices within Elem1 for each face DOF.
                        mfem::Array<int> dofs1;
                        fes_.GetElementDofs(faceTrans->Elem1No, dofs1);
                        std::vector<int> localDofIndices(node_pair_local.first.size());
                        for (size_t p = 0; p < node_pair_local.first.size(); p++) {
                            localDofIndices[p] = dofs1.Find(node_pair_local.first[p]);
                        }

                        auto perDofNormals = computePerDofFaceNormals(
                            *const_cast<mfem::FaceElementTransformations*>(faceTrans),
                            fes_, localDofIndices);

                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            auto dof_normal = perDofNormals[p];
                            if (swap) dof_normal *= -1.0;
                            std::array<double,9> dof_rot{};
                            buildFaceRotationMatrix(dof_normal, mesh->Dimension(), dof_rot.data());

                            if (swap) {
                                raw_pairs[bdr_attr].push_back({
                                    NodePair(node_pair_local.second[p], node_pair_local.first[p]),
                                    dof_rot
                                });
                            } else {
                                raw_pairs[bdr_attr].push_back({
                                    NodePair(node_pair_local.first[p], node_pair_local.second[p]),
                                    dof_rot
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
                        mfem::Array<int> elem_faces, elem_face_ori;
                        if (mesh->Dimension() == 2) {
                            mesh->GetElementEdges(faceTrans->Elem1No, elem_faces, elem_face_ori);
                        } else {
                            mesh->GetElementFaces(faceTrans->Elem1No, elem_faces, elem_face_ori);
                        }
                        const int boundary_face = mesh->GetBdrElementFaceIndex(b);
                        const int local_face = elem_faces.Find(boundary_face);
                        if (local_face == -1) {
                            throw std::runtime_error("Could not map SGBC boundary face to a local element face.");
                        }
                        tagBdrAttributesForSubMesh(local_face, elementSubMesh);
                        mfem::FiniteElementSpace subFES(&elementSubMesh, fec);
                        auto node_pair_local = buildConnectivityForBdrFace(*faceTrans, fes_, subFES);

                        // Compute per-DOF normals on curved faces.
                        mfem::Array<int> dofs1_bdr;
                        fes_.GetElementDofs(faceTrans->Elem1No, dofs1_bdr);
                        std::vector<int> localDofIndices_bdr(node_pair_local.first.size());
                        for (size_t p = 0; p < node_pair_local.first.size(); p++) {
                            localDofIndices_bdr[p] = dofs1_bdr.Find(node_pair_local.first[p]);
                        }

                        auto perDofNormals_bdr = computePerDofFaceNormals(
                            *const_cast<mfem::FaceElementTransformations*>(faceTrans),
                            fes_, localDofIndices_bdr);

                        for (auto p = 0; p < node_pair_local.first.size(); p++){
                            std::array<double,9> dof_rot{};
                            buildFaceRotationMatrix(perDofNormals_bdr[p], mesh->Dimension(), dof_rot.data());
                            raw_pairs[bdr_attr].push_back({
                                NodePair(node_pair_local.first[p], -1),
                                dof_rot
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

                    const ExporterProbe* sgbc_exporter_probe =
                        (sbcps[p].exporter_probe && !probes.exporterProbes.empty())
                            ? &probes.exporterProbes.front()
                            : nullptr;

                    auto wrap = SGBCWrapper::buildSGBCWrapper(sbcps[p], final_time, sgbc_exporter_probe);
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

    Probes emptyProbes; 

    // Keep TFSF submesh infrastructure for planewave source evaluation
    if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {
        srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
        srcmngr_.initDirectPlanewaveEval();
        auto src_sm = static_cast<mfem::SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
        mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), tfsf_sub_to_parent_ids_);
    }

    // Build all operators on the global mesh
    ProblemDescription pd(model_, emptyProbes, srcmngr_.sources, opts_);
    DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, fes_);
    globalOperator_ = dgops.buildGlobalOperator();

    if (model_.hasPML() && model_.getPMLAuxLayout()) {
        PMLOperator_ = dgops.buildPMLOperator(*model_.getPMLAuxLayout());
        initPMLLocalizationDiagnostics();
        if (Mpi::WorldRank() == 0) {
            std::cout << "[PML] Extended ODE size: " << total_state_size_
                      << " (field block " << 6 * fes_.GetNDofs()
                      << " + n_aux " << model_.getPMLAuxLayout()->nAux() << ")"
                      << ", outer_elem=" << pml_outer_element_
                      << std::endl;
        }
        if (isPMLMultProbeEnabled()) {
            runPMLMultProbes();
        }
    }

    if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {
        MPI_Barrier(MPI_COMM_WORLD);
        TFSFOperator_ = dgops.buildTFSFGlobalOperator();
        logTFSFPartitionFaceCoverage(model_.getMesh(),
            model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));

        if (opts_.export_evolution_operator) {
            if (Mpi::WorldSize() == 1) {
                std::filesystem::path export_dir = std::filesystem::path("Exports") / "Operators" / model_.meshName_;
                if (!std::filesystem::exists(export_dir)) {
                    std::filesystem::create_directories(export_dir);
                }
                std::filesystem::path file_path = export_dir / (model_.meshName_ + "_tfsf.csr");
                std::ofstream ofs(file_path);
                if (!ofs.is_open()) {
                    throw std::runtime_error("Could not open file for writing: " + file_path.string());
                }
                TFSFOperator_->PrintCSR2(ofs);
                ofs.close();
                std::cout << "TFSF operator exported to " << file_path << std::endl;

                // Export TFSF mapping matrix T (maps TFSF submesh DOFs to parent mesh DOFs)
                const int tfsf_ndofs = tfsf_sub_to_parent_ids_.Size();
                const int parent_ndofs = fes_.GetNDofs();
                auto mapping_matrix = std::make_unique<mfem::SparseMatrix>(6 * parent_ndofs, 6 * tfsf_ndofs);

                for (int comp = 0; comp < 6; ++comp) {
                    for (int tfsf_dof = 0; tfsf_dof < tfsf_ndofs; ++tfsf_dof) {
                        int parent_dof = tfsf_sub_to_parent_ids_[tfsf_dof];
                        int row = comp * parent_ndofs + parent_dof;
                        int col = comp * tfsf_ndofs + tfsf_dof;
                        mapping_matrix->Set(row, col, 1.0);
                    }
                }
                mapping_matrix->Finalize();

                std::filesystem::path mapping_path = export_dir / (model_.meshName_ + "_tfsf_mapping.csr");
                std::ofstream ofs_map(mapping_path);
                if (!ofs_map.is_open()) {
                    throw std::runtime_error("Could not open file for writing: " + mapping_path.string());
                }
                mapping_matrix->PrintCSR2(ofs_map);
                ofs_map.close();
                std::cout << "TFSF mapping matrix exported to " << mapping_path << std::endl;
            }
        }
    }

    const bool has_rcs_probe = !probes.rcsSurfaceProbes.empty();
    const bool has_mor_probe = !probes.morStateProbes.empty();
    if (opts_.export_evolution_operator && has_rcs_probe && has_mor_probe) {
        if (Mpi::WorldSize() == 1) {
            const auto rcs_tags = collectRCSSurfaceTags(probes);
            if (!rcs_tags.empty()) {
                auto marker = buildSurfaceMarker(rcs_tags, fes_);
                NearToFarFieldSubMesher ntff_submesher(model_.getConstMesh(), fes_, marker);

                auto* ntff_submesh = ntff_submesher.getSubMesh();
                auto dgfec = dynamic_cast<const mfem::DG_FECollection*>(fes_.FEColl());
                if (!dgfec) {
                    throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
                }

                mfem::FiniteElementSpace ntff_fes(ntff_submesh, dgfec);
                mfem::Array<int> farfield_sub_to_parent_ids;
                mfem::SubMeshUtils::BuildVdofToVdofMap(ntff_fes,
                                                       fes_,
                                                       ntff_submesh->GetFrom(),
                                                       ntff_submesh->GetParentElementIDMap(),
                                                       farfield_sub_to_parent_ids);

                const int farfield_ndofs = farfield_sub_to_parent_ids.Size();
                const int parent_ndofs = fes_.GetNDofs();
                auto farfield_mapping_matrix = std::make_unique<mfem::SparseMatrix>(6 * farfield_ndofs, 6 * parent_ndofs);

                for (int comp = 0; comp < 6; ++comp) {
                    for (int sub_dof = 0; sub_dof < farfield_ndofs; ++sub_dof) {
                        const int parent_dof = farfield_sub_to_parent_ids[sub_dof];
                        const int row = comp * farfield_ndofs + sub_dof;
                        const int col = comp * parent_ndofs + parent_dof;
                        farfield_mapping_matrix->Set(row, col, 1.0);
                    }
                }
                farfield_mapping_matrix->Finalize();

                std::filesystem::path export_dir = std::filesystem::path("Exports") / "Operators" / model_.meshName_;
                if (!std::filesystem::exists(export_dir)) {
                    std::filesystem::create_directories(export_dir);
                }

                std::filesystem::path farfield_path = export_dir / (model_.meshName_ + "_farfield.csr");
                std::ofstream ofs_farfield(farfield_path);
                if (!ofs_farfield.is_open()) {
                    throw std::runtime_error("Could not open file for writing: " + farfield_path.string());
                }
                farfield_mapping_matrix->PrintCSR2(ofs_farfield);
                ofs_farfield.close();
                std::cout << "Farfield mapping matrix exported to " << farfield_path << std::endl;
            }
        }
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

static void rotateGlobalFieldsToLocal(
    const double* R,
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields,
    int dof_index, double* out_eh)
{
    double eg[3], hg[3];
    for (int d = 0; d < 3; ++d) {
        eg[d] = fields.get(E, static_cast<Direction>(d))[dof_index];
        hg[d] = fields.get(H, static_cast<Direction>(d))[dof_index];
    }
    for (int d = 0; d < 3; ++d) {
        out_eh[d]     = R[3*d]*eg[0] + R[3*d+1]*eg[1] + R[3*d+2]*eg[2];
        out_eh[3 + d] = R[3*d]*hg[0] + R[3*d+1]*hg[1] + R[3*d+2]*hg[2];
    }
}

GlobalEvolution::~GlobalEvolution() = default;

void GlobalEvolution::commitSGBCCheckpoint(double base_time, double dt,
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields)
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

            // Extract face-local ghost values at time t for interpolation.
            const auto& s = states[i];
            rotateGlobalFieldsToLocal(s.rot, fields, s.global_pair.first, cp[i].ghost_old);
            if (s.global_pair.second != -1) {
                rotateGlobalFieldsToLocal(s.rot, fields, s.global_pair.second, cp[i].ghost_old + 6);
            }
        }
    }
}

void GlobalEvolution::finalizeSGBCStep(
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields)
{
    // Restore from checkpoint and copy interpolation endpoints into working states.
    for (auto& [tag, states] : sgbc_states_) {
        const auto& cp = sgbc_states_checkpoint_.at(tag);
        for (size_t i = 0; i < states.size(); ++i) {
            std::memcpy(states[i].fields_state.GetData(),
                        cp[i].fields_state.GetData(),
                        cp[i].fields_state.Size() * sizeof(double));

            // ghost_old was captured at checkpoint time (t)
            std::memcpy(states[i].ghost_old, cp[i].ghost_old, 12 * sizeof(double));

            // Extract ghost_new from post-RK fields (t + dt)
            const NodePair& pair = states[i].global_pair;
            rotateGlobalFieldsToLocal(states[i].rot, fields, pair.first, states[i].ghost_new);
            if (pair.second != -1) {
                rotateGlobalFieldsToLocal(states[i].rot, fields, pair.second, states[i].ghost_new + 6);
            }
        }
    }

    // Check quiescence: skip entire finalization when all SGBC faces are negligible
    const double thr2 = sgbc_skip_threshold_ * sgbc_skip_threshold_;
    auto ifaceNormSqGF = [&fields](const SGBCState& s) {
        double n2 = 0.0;
        for (int d = 0; d < 3; ++d) {
            double ve = fields.get(E, static_cast<Direction>(d))[s.global_pair.first];
            double vh = fields.get(H, static_cast<Direction>(d))[s.global_pair.first];
            n2 += ve * ve + vh * vh;
            if (s.global_pair.second != -1) {
                ve = fields.get(E, static_cast<Direction>(d))[s.global_pair.second];
                vh = fields.get(H, static_cast<Direction>(d))[s.global_pair.second];
                n2 += ve * ve + vh * vh;
            }
        }
        return n2;
    };

    bool all_quiescent = true;
    for (const auto& [tag, states] : sgbc_states_) {
        for (const auto& state : states) {
            if (state.fields_state.Norml2() >= sgbc_skip_threshold_ ||
                ifaceNormSqGF(state) >= thr2) {
                all_quiescent = false;
                break;
            }
        }
        if (!all_quiescent) break;
    }
    if (all_quiescent) {
        for (auto& [tag, wrapper] : sgbc_wrapper_map_) {
            wrapper->setOldTime(sgbc_step_base_time_ + sgbc_step_dt_);
            wrapper->updateProbes(sgbc_step_base_time_ + sgbc_step_dt_);
        }
        return;
    }

    // Advance SGBC for full dt with linearly interpolated ghost boundary data

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
                double alpha = static_cast<double>(step + 1) / nsteps;
                w->updateFieldsWithInterpolatedGhost(alpha, state);
                w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
            }
            w->saveState(state);
        }
        for (auto& [tag, _] : sgbc_wrapper_map_) {
            _->setOldTime(sgbc_step_base_time_ + dt);
            _->updateProbes(sgbc_step_base_time_ + dt);
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
                    double alpha = static_cast<double>(step + 1) / nsteps;
                    w->updateFieldsWithInterpolatedGhost(alpha, state);
                    w->solve(sgbc_step_base_time_ + step * actual_sub_dt, actual_sub_dt);
                }
                w->saveState(state);
            }
            w->setOldTime(sgbc_step_base_time_ + dt);
            w->updateProbes(sgbc_step_base_time_ + dt);
        }
    }
}

void GlobalEvolution::applyTFSFSourceToVector(double t_stage, int ndofs, int nbrDofs,
                                              mfem::Vector& result_vector) const
{
    if (!TFSFOperator_) return;
    if (tfsfSourceIndices_.empty()) return;

    auto *tfsf_space = srcmngr_.getGlobalTFSFSpace();
    if (!tfsf_space) return;

    const int blockSize = ndofs + nbrDofs;

    // Ensure global-layout planewave vector is allocated (first call only)
    if (multWorkVec_.Size() != 6 * blockSize) {
        multWorkVec_.SetSize(6 * blockSize);
        multWorkVec_.UseDevice(true);
        multWorkVec_ = 0.0;
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

    // Check if evaluated planewave is negligible — skip scatter + SpMV.
    // Not permanent: CW sources (e.g. cosine) pass through zero between lobes,
    // so we re-evaluate every call and only skip the current one.
    {
        double norm2 = 0.0;
        for (int f : {E, H}) {
            for (int d : {X, Y, Z}) {
                double n = func[f][d].Norml2();
                norm2 += n * n;
            }
        }
        if (norm2 < tfsf_skip_threshold_ * tfsf_skip_threshold_) {
            return;
        }
    }

    // Scatter submesh planewave values into multWorkVec_ (local DOFs only;
    // partitioner guarantees all TFSF faces are rank-local, no ghost exchange needed).
#ifdef SEMBA_DGTD_ENABLE_CUDA
    scatter_tfsf_to_assembled_gpu(tfsf_sub_to_parent_ids_, func, multWorkVec_, blockSize);
#else
    for (int v = 0; v < tfsf_sub_to_parent_ids_.Size(); ++v) {
        const int parent_dof = tfsf_sub_to_parent_ids_[v];
        for (int d : {X, Y, Z}) {
            multWorkVec_[d * blockSize + parent_dof] = func[E][d][v];
            multWorkVec_[(3 + d) * blockSize + parent_dof] = func[H][d][v];
        }
    }
#endif

    // Apply TFSF operator and subtract from result: out -= TFSFOp * planewave
    TFSFOperator_->AddMult(multWorkVec_, result_vector, -1.0);

    // Restore zeroed state for next call
#ifdef SEMBA_DGTD_ENABLE_CUDA
    zero_tfsf_assembled_gpu(tfsf_sub_to_parent_ids_, multWorkVec_, blockSize);
#else
    for (int d = X; d <= Z; ++d) {
        std::memset(multWorkVec_.GetData() + d * blockSize, 0, blockSize * sizeof(double));
        std::memset(multWorkVec_.GetData() + (3 + d) * blockSize, 0, blockSize * sizeof(double));
    }
#endif
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

    // Ensure ODE state is readable (device→host if needed).
    (void)in.Read();

    // 1) SGBC sub-solve FIRST — independent of neighbor data, overlaps with
    //    other ranks' local work so they don't idle-wait during exchange.
#ifdef SHOW_TIMER_INFORMATION
    timerSGBC.Start();
#endif
    bool sgbc_all_quiescent = false;
    if (SGBCOperator_ && !sgbcWrappers_.empty()){

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

        // Helper: compute squared norm of global interface DOFs for a given state.
        auto ifaceNormSq = [&in, ndofs](const SGBCState& s) {
            double n2 = 0.0;
            for (int d = 0; d < 6; ++d) {
                double v = in[d * ndofs + s.global_pair.first];
                n2 += v * v;
                if (s.global_pair.second != -1) {
                    v = in[d * ndofs + s.global_pair.second];
                    n2 += v * v;
                }
            }
            return n2;
        };
        const double thr2 = sgbc_skip_threshold_ * sgbc_skip_threshold_;

        // Check if ALL SGBC faces are quiescent (skip sub-solve + flux injection)
        sgbc_all_quiescent = true;
        for (const auto& [tag, states] : sgbc_states_) {
            for (const auto& state : states) {
                if (state.fields_state.Norml2() >= sgbc_skip_threshold_ ||
                    ifaceNormSq(state) >= thr2) {
                    sgbc_all_quiescent = false;
                    break;
                }
            }
            if (!sgbc_all_quiescent) break;
        }

        // Advance SGBC to stage time using intermediate global fields as ghost data
        if (!sgbc_all_quiescent && dt_to_stage > 1e-12) {
#ifdef SEMBA_DGTD_ENABLE_OPENMP
            if (!sgbc_thread_pool_.empty()) {
                #pragma omp parallel for schedule(dynamic) num_threads(sgbc_omp_threads_)
                for (size_t ti = 0; ti < sgbc_tasks_.size(); ++ti) {
                    const auto& task = sgbc_tasks_[ti];
                    auto& state = sgbc_states_[task.tag][task.state_index];

                    // Skip sub-solve when fields are negligible
                    double stateN2 = 0.0;
                    for (int j = 0; j < state.fields_state.Size(); ++j) {
                        stateN2 += state.fields_state[j] * state.fields_state[j];
                    }
                    if (stateN2 < thr2 && ifaceNormSq(state) < thr2) continue;

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
                        // Skip sub-solve when fields are negligible
                        if (state.fields_state.Norml2() < sgbc_skip_threshold_ &&
                            ifaceNormSq(state) < thr2) continue;

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
    {
        const double* in_data = in.Read();
        for (int d = X; d <= Z; ++d) {
            double* e_data = eOld_[d].Write();
            double* h_data = hOld_[d].Write();
            std::memcpy(e_data, in_data + d * ndofs, ndofs * sizeof(double));
            std::memcpy(h_data, in_data + (3 + d) * ndofs, ndofs * sizeof(double));
        }
    }
#endif
    for (int d = X; d <= Z; ++d) {
        eOld_[d].ExchangeFaceNbrData();
        hOld_[d].ExchangeFaceNbrData();
    }
#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (mfem::Device::Allows(mfem::Backend::CUDA) && nbrDofs > 0) {
        sync_cuda_face_nbr_halos(eOld_, hOld_);
    }
#endif
#ifdef SHOW_TIMER_INFORMATION
    timerExchange.Stop();
#endif

    // 3) Assemble multWorkVec_ (local + neighbor data)
    if (multWorkVec_.Size() != 6 * blockSize) {
        multWorkVec_.SetSize(6 * blockSize);
        multWorkVec_.UseDevice(true);
    }

#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (mfem::Device::Allows(mfem::Backend::CUDA)) {
        load_eh_to_innew_gpu(in, multWorkVec_, ndofs, nbrDofs);
        if (nbrDofs > 0) {
            load_nbr_to_innew_gpu(eOld_, hOld_, multWorkVec_, ndofs, nbrDofs);
        }
    } else
#endif
    {
        for (int d = X; d <= Z; ++d) {
            multWorkVec_.SetVector(eOld_[d],       d      * blockSize);
            multWorkVec_.SetVector(hOld_[d],  (3 + d)     * blockSize);
        }
        if (nbrDofs) {
            for (int d = X; d <= Z; ++d) {
                mfem::Vector &eNbr = eOld_[d].FaceNbrData();
                mfem::Vector &hNbr = hOld_[d].FaceNbrData();
                multWorkVec_.SetVector(eNbr,      d      * blockSize + ndofs);
                multWorkVec_.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
            }
        }
    }

    // 4) Apply globalOperator_ (DG flux) 
#ifdef SHOW_TIMER_INFORMATION
    timerApplyA.Start();
#endif
    if (out.Size() != total_state_size_) {
        out.SetSize(total_state_size_);
        out.UseDevice(true);
    }
    out = 0.0;
    if (!shouldSkipGlobalOperatorInMult()) {
        globalOperator_->Mult(multWorkVec_, out);
    }
#ifdef SHOW_TIMER_INFORMATION
    timerApplyA.Stop();
#endif

    if (!shouldSkipPMLOperatorInMult()) {
        applyPMLCoupling(in, out);
    }

    // S3: Zero multWorkVec_ so it can be reused for TFSF and SGBC injection.
    multWorkVec_ = 0.0;

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
    //    Skip entirely when all SGBC faces are quiescent.
    if (SGBCOperator_ && !sgbcWrappers_.empty() && !sgbc_all_quiescent){
        // multWorkVec_ is already sized and zeroed (after step 4 zeroing + TFSF cleanup).
        // On GPU builds 'in' is device-authoritative after steps 2-5; HostRead()
        // syncs device→host once for all operator[]-based DOF reads below.
        // multWorkVec_ was zeroed on device; HostReadWrite() syncs it to host so that
        // operator[] assignments write to the correct (host) buffer.
        in.HostRead();
        multWorkVec_.HostReadWrite();

        // Helper: fill multWorkVec_ slot for a given BC type.
        // For PEC/PMC/SMA: mirror/absorb from global field at source_dof.
        // For transparent (SGBC): rotate slab interior field at slab_idx to global frame.
        auto fillBCSlot = [&](BdrCond bdr, int source_dof, int target_slot,
                              const double* rot, const double* fields_state,
                              int local_size, int slab_idx) {
            if (bdr == BdrCond::PEC) {
                for (int d = X; d <= Z; ++d) {
                    multWorkVec_[d * blockSize + target_slot]       = -in[d * ndofs + source_dof];
                    multWorkVec_[(3 + d) * blockSize + target_slot] =  in[(3 + d) * ndofs + source_dof];
                }
            } else if (bdr == BdrCond::PMC) {
                for (int d = X; d <= Z; ++d) {
                    multWorkVec_[d * blockSize + target_slot]       =  in[d * ndofs + source_dof];
                    multWorkVec_[(3 + d) * blockSize + target_slot] = -in[(3 + d) * ndofs + source_dof];
                }
            } else if (bdr == BdrCond::SMA) {
                for (int d = X; d <= Z; ++d) {
                    multWorkVec_[d * blockSize + target_slot]       = in[d * ndofs + source_dof];
                    multWorkVec_[(3 + d) * blockSize + target_slot] = in[(3 + d) * ndofs + source_dof];
                }
            } else {
                // Transparent: rotate slab local→global
                double el[3], hl[3];
                for (int d = 0; d < 3; ++d) {
                    el[d] = fields_state[d * local_size + slab_idx];
                    hl[d] = fields_state[(3 + d) * local_size + slab_idx];
                }
                for (int d = 0; d < 3; ++d) {
                    multWorkVec_[d * blockSize + target_slot]       = rot[d]*el[0] + rot[3+d]*el[1] + rot[6+d]*el[2];
                    multWorkVec_[(3 + d) * blockSize + target_slot] = rot[d]*hl[0] + rot[3+d]*hl[1] + rot[6+d]*hl[2];
                }
            }
        };

        // Helper: copy global field values into multWorkVec_ slot
        auto copyGlobalSlot = [&](int dof, int slot) {
            for (int d = X; d <= Z; ++d) {
                multWorkVec_[d * blockSize + slot] = in[d * ndofs + dof];
                multWorkVec_[(3 + d) * blockSize + slot] = in[(3 + d) * ndofs + dof];
            }
        };

        for (auto& [tag, states] : sgbc_states_) {
            auto it = sgbc_wrapper_map_.find(tag);
            if (it == sgbc_wrapper_map_.end()) continue;
            auto w = it->second;
            int local_size = w->getLocalFieldSize();
            int idx_left = w->getLeftInterfaceIndex();
            const auto& left_info = w->getProperties().sgbc_bdr_info.first;
            const auto left_bdr = left_info.isOn ? left_info.bdrCond : BdrCond::SGBC;

            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                if (gr == -1) {
                    // Boundary SGBC: the source operator acts on the
                    // SGBC trace directly at the boundary DOF slot.
                    fillBCSlot(left_bdr, gl, gl, state.rot,
                               state.fields_state.GetData(), local_size, idx_left);
                } else {
                    // Elem1 slot: global field value (u_L) — already in global frame
                    copyGlobalSlot(gl, gl);
                    // Elem2 slot: what Elem1 sees from the slab's left side
                    fillBCSlot(left_bdr, gl, gr, state.rot,
                               state.fields_state.GetData(), local_size, idx_left);
                }
            }
        }
        // Targeted SpMV: only compute rows for Elem1 DOFs, add directly to out.
        {
            const double* sv = multWorkVec_.HostRead();
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
                    multWorkVec_[d * blockSize + gl] = 0.0;
                    multWorkVec_[(3 + d) * blockSize + gl] = 0.0;
                    if (gr != -1) {
                        multWorkVec_[d * blockSize + gr] = 0.0;
                        multWorkVec_[(3 + d) * blockSize + gr] = 0.0;
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
            const auto& right_info = w->getProperties().sgbc_bdr_info.second;
            const auto right_bdr = right_info.isOn ? right_info.bdrCond : BdrCond::SGBC;

            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                if (gr == -1) continue;

                fillBCSlot(right_bdr, gr, gl, state.rot,
                           state.fields_state.GetData(), local_size, idx_right);
                // Elem2 slot: global field value (u_R) — already in global frame
                copyGlobalSlot(gr, gr);
            }
        }
        // Targeted SpMV: only compute rows for Elem2 DOFs, add directly to out.
        {
            const double* sv = multWorkVec_.HostRead();
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

        // Clear entries written in Pass 2 so multWorkVec_ stays zeroed for next call
        for (auto& [tag, states] : sgbc_states_) {
            for (const auto& state : states) {
                int gl = state.global_pair.first;
                int gr = state.global_pair.second;
                if (gr == -1) continue;
                for (int d = X; d <= Z; ++d) {
                    multWorkVec_[d * blockSize + gl] = 0.0;
                    multWorkVec_[(3 + d) * blockSize + gl] = 0.0;
                    multWorkVec_[d * blockSize + gr] = 0.0;
                    multWorkVec_[(3 + d) * blockSize + gr] = 0.0;
                }
            }
        }
#ifdef SEMBA_DGTD_ENABLE_CUDA
        // Steps 3-5 wrote 'out' on device; the SGBC SpMV passes above wrote to
        // the host copy via HostReadWrite(). Sync host→device so the ODE
        // integrator sees a consistent device-authoritative 'out'.
        (void)out.ReadWrite();
#endif
    }

    // PML (and possibly TFSF/SGBC host writes) may have updated the host copy of
    // 'out'; sync to device so the RK stage vectors see a consistent rate vector.
    (void)out.ReadWrite();

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

void GlobalEvolution::runPMLMultProbes()
{
    if (!PMLOperator_ || Mpi::WorldRank() != 0) {
        return;
    }

    const int ndofs = fes_.GetNDofs();
    const PMLAuxLayout* layout = model_.getPMLAuxLayout();

    std::unordered_set<int> pml_dofs;
    for (int d : pml_inner_pml_dof_indices_) {
        pml_dofs.insert(d);
    }
    for (int d : pml_outer_dof_indices_) {
        pml_dofs.insert(d);
    }

    int vacuum_ey_dof = -1;
    for (int i = 0; i < ndofs; ++i) {
        if (pml_dofs.find(i) == pml_dofs.end()) {
            vacuum_ey_dof = i;
            break;
        }
    }
    const int pml_ey_dof =
        pml_inner_pml_dof_indices_.empty() ? -1 : pml_inner_pml_dof_indices_.front();

    int psi_index = -1;
    if (layout && layout->numStretchDirections() > 0) {
        const Direction stretch = layout->stretchDirection(0);
        const int psi_off = layout->psiEOffset(stretch, Z);
        const int psi_dof = (pml_ey_dof >= 0) ? pml_ey_dof : 0;
        psi_index = psi_off + psi_dof;
    }

    mfem::Vector in(total_state_size_);
    mfem::Vector out(total_state_size_);
    const int hz_base = 3 * ndofs + static_cast<int>(Z);
    const int ey_base = ndofs * static_cast<int>(Y);
    const int nbrDofs = fes_.num_face_nbr_dofs;
    const int blockSize = ndofs + nbrDofs;
    if (multWorkVec_.Size() != 6 * blockSize) {
        multWorkVec_.SetSize(6 * blockSize);
    }

    auto run_probe = [&](const char* name, int ey_dof, int psi_row) {
        in = 0.0;
        out = 0.0;
        multWorkVec_ = 0.0;
        if (ey_dof >= 0) {
            in[ey_base + ey_dof] = 1.0;
            multWorkVec_[ey_base + ey_dof] = 1.0;
        }
        if (psi_row >= 0) {
            in[psi_row] = 1.0;
        }
        if (total_state_size_ > 6 * ndofs) {
            std::memset(out.GetData() + 6 * ndofs, 0,
                        static_cast<size_t>(total_state_size_ - 6 * ndofs) * sizeof(double));
        }

        double max_hz_global = 0.0;
        double max_hz_after_pml = 0.0;
        double max_psi_dot = 0.0;

        globalOperator_->Mult(multWorkVec_, out);
        for (int i = 0; i < ndofs; ++i) {
            max_hz_global = std::max(max_hz_global, std::abs(out[hz_base + i]));
        }

        applyPMLCoupling(in, out);
        for (int i = 0; i < ndofs; ++i) {
            max_hz_after_pml = std::max(max_hz_after_pml, std::abs(out[hz_base + i]));
        }
        for (int i = 6 * ndofs; i < total_state_size_; ++i) {
            max_psi_dot = std::max(max_psi_dot, std::abs(out[i]));
        }
        pmlAuditLog(std::string(name) + ": max|dHz/dt| global=" +
                    std::to_string(max_hz_global) + " after_pml=" +
                    std::to_string(max_hz_after_pml) + " max|dpsi/dt|=" +
                    std::to_string(max_psi_dot));
    };

    pmlAuditLog("=== PML Mult probes (globalOperator + applyPMLCoupling only) ===");
    if (vacuum_ey_dof >= 0) {
        run_probe("P1_vacuum_Ey", vacuum_ey_dof, -1);
    }
    if (pml_ey_dof >= 0) {
        run_probe("P2_pml_Ey", pml_ey_dof, -1);
    }
    if (psi_index >= 0) {
        run_probe("P3_unit_psiEz", -1, psi_index);
    }
    run_probe("P4_zero", -1, -1);
}

void GlobalEvolution::initPMLLocalizationDiagnostics()
{
    pml_outer_element_ = -1;
    pml_outer_dof_indices_.clear();
    pml_inner_pml_dof_indices_.clear();
    if (!PMLOperator_) {
        return;
    }

    mfem::Mesh& mesh = *fes_.GetMesh();
    const mfem::Array<int> pml_marker = model_.buildPMLVolumeMarker();
    const int dim = mesh.Dimension();
    double best_coord = -1e300;
    std::vector<int> pml_elements;

    for (int el = 0; el < mesh.GetNE(); ++el) {
        const int attr = mesh.GetAttribute(el);
        if (attr < 1 || attr > pml_marker.Size() || pml_marker[attr - 1] == 0) {
            continue;
        }
        pml_elements.push_back(el);
        const mfem::Geometry::Type geom = mesh.GetElementBaseGeometry(el);
        const mfem::IntegrationPoint& ip = mfem::Geometries.GetCenter(geom);
        mfem::ElementTransformation* tr = mesh.GetElementTransformation(el);
        tr->SetIntPoint(&ip);
        mfem::Vector x(dim);
        tr->Transform(ip, x);
        const double coord = x(0);
        if (coord > best_coord) {
            best_coord = coord;
            pml_outer_element_ = el;
        }
    }

    const mfem::Table& elem_dof_table = fes_.GetElementToDofTable();
    auto collect_element_dofs = [&](int el, std::unordered_set<int>& dof_set) {
        const int* dofs = elem_dof_table.GetRow(el);
        const int n = elem_dof_table.RowSize(el);
        for (int i = 0; i < n; ++i) {
            dof_set.insert(dofs[i]);
        }
    };

    std::unordered_set<int> outer_set;
    std::unordered_set<int> inner_set;
    if (pml_outer_element_ >= 0) {
        collect_element_dofs(pml_outer_element_, outer_set);
    }
    for (int el : pml_elements) {
        if (el == pml_outer_element_) {
            continue;
        }
        collect_element_dofs(el, inner_set);
    }
    for (int dof : outer_set) {
        pml_outer_dof_indices_.push_back(dof);
    }
    for (int dof : inner_set) {
        if (outer_set.find(dof) == outer_set.end()) {
            pml_inner_pml_dof_indices_.push_back(dof);
        }
    }
    std::sort(pml_outer_dof_indices_.begin(), pml_outer_dof_indices_.end());
    std::sort(pml_inner_pml_dof_indices_.begin(), pml_inner_pml_dof_indices_.end());
}

void GlobalEvolution::applyPMLCoupling(const mfem::Vector& in, mfem::Vector& out) const
{
    if (!PMLOperator_) {
        return;
    }

    const int ndofs = fes_.GetNDofs();
    const int nbrDofs = fes_.num_face_nbr_dofs;
    const int blockSize = ndofs + nbrDofs;
    const int n_aux = total_state_size_ - 6 * ndofs;
    const int work_size = 6 * blockSize + n_aux;

    if (pmlWorkVec_.Size() != work_size) {
        pmlWorkVec_.SetSize(work_size);
        pmlWorkVec_.UseDevice(true);
    }

    std::memcpy(
        pmlWorkVec_.HostWrite(),
        multWorkVec_.HostRead(),
        static_cast<size_t>(6 * blockSize) * sizeof(double));
    std::memcpy(
        pmlWorkVec_.HostWrite() + 6 * blockSize,
        in.HostRead() + 6 * ndofs,
        static_cast<size_t>(n_aux) * sizeof(double));

    static bool pml_workvec_checked = false;
    if (isPMLOperatorAuditEnabled() && !pml_workvec_checked && Mpi::WorldRank() == 0) {
        pml_workvec_checked = true;
        const double* mw = multWorkVec_.HostRead();
        const double* pw = pmlWorkVec_.HostRead();
        double max_diff = 0.0;
        for (int i = 0; i < 6 * blockSize; ++i) {
            max_diff = std::max(max_diff, std::abs(mw[i] - pw[i]));
        }
        pmlAuditLog("multWorkVec vs pmlWorkVec fields max|diff|=" +
                    std::to_string(max_diff));
    }

    const double* out_before = out.HostRead();
    double max_hz_before = 0.0;
    const int hz_base = 3 * ndofs + static_cast<int>(Z);
    if (isPMLOperatorAuditEnabled()) {
        for (int dof : pml_outer_dof_indices_) {
            max_hz_before = std::max(max_hz_before, std::abs(out_before[hz_base + dof]));
        }
    }

    PMLOperator_->AddMult(pmlWorkVec_, out, getPMLOperatorMultSign());

    if (isPMLOperatorAuditEnabled() && Mpi::WorldRank() == 0 &&
        pml_diag_mult_count_ == 1 && !pml_outer_dof_indices_.empty()) {
        const double* out_after = out.HostRead();
        double max_hz_after = 0.0;
        for (int dof : pml_outer_dof_indices_) {
            max_hz_after = std::max(max_hz_after, std::abs(out_after[hz_base + dof]));
        }
        pmlAuditLog("outer PML max|dHz/dt| before_pml=" + std::to_string(max_hz_before) +
                    " after_pml_add=" + std::to_string(max_hz_after));
    }

    ++pml_diag_mult_count_;
    if (Mpi::WorldRank() == 0 &&
        (pml_diag_mult_count_ == 1 || pml_diag_mult_count_ % 500 == 0)) {
        const double* in_data = in.HostRead();
        const double* out_data = out.HostRead();
        const PMLAuxLayout* layout = model_.getPMLAuxLayout();
        double max_psi = 0.0;
        double max_psi_dot = 0.0;
        double max_psi_outer = 0.0;
        double max_psi_inner = 0.0;
        double max_e_corr = 0.0;
        double max_h_corr = 0.0;
        double max_ey = 0.0;
        double max_ey_outer = 0.0;
        const int ey_base = ndofs * static_cast<int>(Y);

        auto max_psi_on_dofs = [&](int offset, const std::vector<int>& dofs, double& target) {
            for (int dof : dofs) {
                const double v = std::abs(in_data[offset + dof]);
                target = std::max(target, v);
            }
        };

        for (int slot = 0; slot < layout->numStretchDirections(); ++slot) {
            const Direction d = layout->stretchDirection(slot);
            for (Direction c = X; c <= Z; ++c) {
                if (!pmlPsiEComponentActive(c, d)) {
                    continue;
                }
                const int off_e = layout->psiEOffset(d, c);
                for (int i = 0; i < ndofs; ++i) {
                    max_psi = std::max(max_psi, std::abs(in_data[off_e + i]));
                    max_psi_dot = std::max(max_psi_dot, std::abs(out_data[off_e + i]));
                }
                max_psi_on_dofs(off_e, pml_outer_dof_indices_, max_psi_outer);
                max_psi_on_dofs(off_e, pml_inner_pml_dof_indices_, max_psi_inner);
            }
            for (Direction c = X; c <= Z; ++c) {
                if (!pmlPsiHComponentActive(c, d)) {
                    continue;
                }
                const int off_h = layout->psiHOffset(d, c);
                for (int i = 0; i < ndofs; ++i) {
                    max_psi = std::max(max_psi, std::abs(in_data[off_h + i]));
                    max_psi_dot = std::max(max_psi_dot, std::abs(out_data[off_h + i]));
                }
                max_psi_on_dofs(off_h, pml_outer_dof_indices_, max_psi_outer);
                max_psi_on_dofs(off_h, pml_inner_pml_dof_indices_, max_psi_inner);
            }
        }
        for (int i = 0; i < 3 * ndofs; ++i) {
            max_e_corr = std::max(max_e_corr, std::abs(out_data[i]));
            max_h_corr = std::max(max_h_corr, std::abs(out_data[3 * ndofs + i]));
        }
        for (int i = 0; i < ndofs; ++i) {
            const double ey = std::abs(in_data[ey_base + i]);
            max_ey = std::max(max_ey, ey);
        }
        for (int dof : pml_outer_dof_indices_) {
            max_ey_outer = std::max(max_ey_outer, std::abs(in_data[ey_base + dof]));
        }
        std::cout << "[PML Mult diag @ call " << pml_diag_mult_count_ << "] "
                  << "max|psi|=" << max_psi
                  << " outer=" << max_psi_outer
                  << " inner_pml=" << max_psi_inner
                  << " max|dpsi/dt|=" << max_psi_dot
                  << " max|out E|=" << max_e_corr
                  << " max|out H|=" << max_h_corr
                  << " max|Ey|=" << max_ey
                  << " outer_Ey=" << max_ey_outer
                  << std::endl;
    }
}

void GlobalEvolution::ImplicitSolve(const double dt,
                                    const mfem::Vector& x,
                                    mfem::Vector& k)
{
    const int n         = x.Size();
    const int ndofs     = fes_.GetNDofs();
    const int nbrDofs   = fes_.num_face_nbr_dofs;
    const int blockSize = ndofs + nbrDofs;
    if (total_state_size_ > 6 * ndofs) {
        MFEM_ABORT("PML extended state is not supported with implicit integrators in Milestone A.");
    }
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
#ifdef SEMBA_DGTD_ENABLE_CUDA
        if (mfem::Device::Allows(mfem::Backend::CUDA)) {
            load_in_to_eh_gpu(x, eOld_, hOld_, ndofs);
        } else
#endif
        {
            const double* x_host = x.HostRead();
            for (int d = X; d <= Z; ++d) {
                double* e_data = eOld_[d].HostWrite();
                double* h_data = hOld_[d].HostWrite();
                std::memcpy(e_data, x_host + d * ndofs, ndofs * sizeof(double));
                std::memcpy(h_data, x_host + (3 + d) * ndofs, ndofs * sizeof(double));
            }
        }
        for (int d = X; d <= Z; ++d) {
            eOld_[d].ExchangeFaceNbrData();
            hOld_[d].ExchangeFaceNbrData();
        }
#ifdef SEMBA_DGTD_ENABLE_CUDA
        if (mfem::Device::Allows(mfem::Backend::CUDA)) {
            if (nbrDofs > 0) {
                sync_cuda_face_nbr_halos(eOld_, hOld_);
            }
            load_eh_to_innew_gpu(x, implicit_inNew_, ndofs, nbrDofs);
            if (nbrDofs > 0) {
                load_nbr_to_innew_gpu(eOld_, hOld_, implicit_inNew_, ndofs, nbrDofs);
            }
        } else
#endif
        {
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
        }
        globalOperator_->Mult(implicit_inNew_, implicit_rhs_);
    }

    // Add TFSF source contribution (skip if no active TFSF sources)
    if (TFSFOperator_ && !tfsfSourceIndices_.empty()) {
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
#ifdef SEMBA_DGTD_ENABLE_CUDA
            if (mfem::Device::Allows(mfem::Backend::CUDA)) {
                load_in_to_eh_gpu(u, eOld_, hOld_, ndofs);
            } else
#endif
            {
                const double* u_host = u.HostRead();
                for (int d = X; d <= Z; ++d) {
                    double* e_data = eOld_[d].HostWrite();
                    double* h_data = hOld_[d].HostWrite();
                    std::memcpy(e_data, u_host + d * ndofs, ndofs * sizeof(double));
                    std::memcpy(h_data, u_host + (3 + d) * ndofs, ndofs * sizeof(double));
                }
            }
            for (int d = X; d <= Z; ++d) {
                eOld_[d].ExchangeFaceNbrData();
                hOld_[d].ExchangeFaceNbrData();
            }
#ifdef SEMBA_DGTD_ENABLE_CUDA
            if (mfem::Device::Allows(mfem::Backend::CUDA)) {
                if (nbrDofs > 0) {
                    sync_cuda_face_nbr_halos(eOld_, hOld_);
                }
                load_eh_to_innew_gpu(u, implicit_inNew_, ndofs, nbrDofs);
                if (nbrDofs > 0) {
                    load_nbr_to_innew_gpu(eOld_, hOld_, implicit_inNew_, ndofs, nbrDofs);
                }
            } else
#endif
            {
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