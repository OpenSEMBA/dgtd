#include "SourcesManager.h"
#include <unordered_map>
#include <unordered_set>

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(const Sources& srcs, mfem::ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>& fields) :
    sources{ srcs }, 
    fes_{ fes }
{
    setInitialFields(fields);
}

void SourcesManager::setInitialFields(Fields<ParFiniteElementSpace, ParGridFunction>& fields)
{
    for (const auto& source : sources) {
        auto src{ dynamic_cast<InitialField*>(source.get()) };
        if (src == nullptr) {
            continue;
        }
        for (FieldType ft: {E, H}) {
            for (auto x : { X, Y, Z }) {
                auto f = [src, ft, x](const Source::Position& pos, Source::Time t) {
                    return src->eval(pos, t, ft, x);
                };
                FunctionCoefficient fc(f);
                GridFunction gf(fields.get(ft, x).FESpace());
                gf.UseDevice(true);
                gf.ProjectCoefficient(fc);
                fields.get(ft, x) += gf;
            }
        }
    }
}

FieldGridFuncs SourcesManager::evalTimeVarField(const Time time, FiniteElementSpace* fes)
{
    if (!cached_tfsf_fields_init_) {
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                cached_tfsf_fields_[ft][d].UseDevice(true);
                cached_tfsf_fields_[ft][d].SetSpace(fes);
            }
        }
        cached_tfsf_fields_init_ = true;
    }

    for (auto ft : { E, H }) {
        for (auto d : { X, Y, Z }) {
            cached_tfsf_fields_[ft][d] = 0.0;
        }
    }

    for (const auto& source : sources) {
        auto tf = dynamic_cast<TotalField*>(source.get());
        if (tf == nullptr) {
            continue;
        }
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                auto f = [tf, ft, d](const Source::Position& pos, Source::Time t) {
                    return tf->eval(pos, t, ft, d);
                };
                FunctionCoefficient func(f);
                func.SetTime(time);
                cached_tfsf_fields_[ft][d].ProjectCoefficient(func);
            }
        }
    }
    return cached_tfsf_fields_;
}

void SourcesManager::markDoFSforTForSF(FieldGridFuncs& gfs, bool isTF)
{
    auto global_tfsf_map = tfsf_submesher_.getGlobalTFSFSubMesh()->GetParentElementIDMap();
    Array<int> secondary_map;
    switch (isTF) {
    case true:
        if (tfsf_submesher_.getSFSubMesh() != NULL) {
            secondary_map = tfsf_submesher_.getSFSubMesh()->GetParentElementIDMap();
        }
        break;
    case false:
        if (tfsf_submesher_.getTFSubMesh() != NULL) {
            secondary_map = tfsf_submesher_.getTFSubMesh()->GetParentElementIDMap();
        }
        break;
    }

    for (int e = 0; e < secondary_map.Size(); e++) {
        Array<int> dofs;
        global_tfsf_fes_->GetElementDofs(global_tfsf_map.Find(secondary_map[e]), dofs);
        for (int i = 0; i < dofs.Size(); i++) {
            for (auto f : { E, H }) {
                for (auto d{ X }; d <= Z; d++) {
                    gfs[f][d].UseDevice(true);
                    gfs[f][d][dofs[i]] = 0.0;
                }
            }
        }
    }
}

void SourcesManager::initTFSFPreReqs(const ParMesh& m, const Array<int>& marker)
{
    initTFSFSubMesher(m, marker);
    initTFSFSpaces();
}

void SourcesManager::initTFSFSubMesher(const ParMesh& m, const Array<int>& marker)
{
    auto sm = TotalFieldScatteredFieldSubMesher(m, marker);
    tfsf_submesher_ = std::move(sm);
}

void SourcesManager::initTFSFSpaces()
{
    if (tfsf_submesher_.getTFSubMesh() != NULL) {
        tf_fes_ = std::make_unique<FiniteElementSpace>(tfsf_submesher_.getTFSubMesh(), fes_.FEColl());
    }
    if (tfsf_submesher_.getSFSubMesh() != NULL) {
        sf_fes_ = std::make_unique<FiniteElementSpace>(tfsf_submesher_.getSFSubMesh(), fes_.FEColl());
    }
    global_tfsf_fes_ = std::make_unique<FiniteElementSpace>(tfsf_submesher_.getGlobalTFSFSubMesh(), fes_.FEColl());
}

void SourcesManager::initDirectPlanewaveEval()
{
    if (!global_tfsf_fes_) return;
    auto* fes = global_tfsf_fes_.get();
    const int ndofs = fes->GetNDofs();
    auto* mesh = fes->GetMesh();

    // 1) Pre-compute physical coordinates for each DOF
    dof_coords_.resize(ndofs);
    for (int el = 0; el < fes->GetNE(); ++el) {
        const auto* fe = fes->GetFE(el);
        auto* T = fes->GetElementTransformation(el);
        const auto& ir = fe->GetNodes();
        mfem::Array<int> dofs;
        fes->GetElementDofs(el, dofs);
        for (int i = 0; i < dofs.Size(); ++i) {
            const auto& ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);
            mfem::Vector coords;
            T->Transform(ip, coords);
            // Pad to 3D for planewave dot-product
            dof_coords_[dofs[i]].SetSize(3);
            dof_coords_[dofs[i]] = 0.0;
            for (int c = 0; c < coords.Size(); ++c) {
                dof_coords_[dofs[i]][c] = coords[c];
            }
        }
    }

    // 2) Pre-compute TF/SF sign mask
    tfsf_sign_.assign(ndofs, 1.0);  // Default: no SF submesh, scale = 1.0
    if (tfsf_submesher_.getSFSubMesh() != NULL) {
        // Build set of SF DOF indices (on the global TFSF submesh)
        auto global_tfsf_map = tfsf_submesher_.getGlobalTFSFSubMesh()->GetParentElementIDMap();

        // Build reverse lookup: parent_elem_id -> global_tfsf_element_index
        std::unordered_map<int, int> parent_to_global;
        for (int i = 0; i < global_tfsf_map.Size(); ++i) {
            parent_to_global[global_tfsf_map[i]] = i;
        }

        // All DOFs start as TF (+0.5)
        tfsf_sign_.assign(ndofs, 0.5);

        // Mark SF DOFs as -0.5
        auto sf_map = tfsf_submesher_.getSFSubMesh()->GetParentElementIDMap();
        std::unordered_set<int> sf_dof_set;
        for (int e = 0; e < sf_map.Size(); ++e) {
            auto it = parent_to_global.find(sf_map[e]);
            if (it == parent_to_global.end()) continue;
            mfem::Array<int> dofs;
            fes->GetElementDofs(it->second, dofs);
            for (int i = 0; i < dofs.Size(); ++i) {
                sf_dof_set.insert(dofs[i]);
            }
        }
        for (int dof : sf_dof_set) {
            tfsf_sign_[dof] = -0.5;
        }
    }

    // 3) Initialize cached grid functions if not already done
    if (!cached_tfsf_fields_init_) {
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                cached_tfsf_fields_[ft][d].UseDevice(true);
                cached_tfsf_fields_[ft][d].SetSpace(fes);
            }
        }
        cached_tfsf_fields_init_ = true;
    }

    direct_eval_ready_ = true;
}

void SourcesManager::evalTimeVarFieldDirect(Time time)
{
    auto* fes = global_tfsf_fes_.get();
    const int ndofs = fes->GetNDofs();

    // Zero all fields
    for (auto ft : { E, H }) {
        for (auto d : { X, Y, Z }) {
            cached_tfsf_fields_[ft][d] = 0.0;
        }
    }

    // Evaluate each TotalField source directly at precomputed DOF positions
    for (const auto& source : sources) {
        auto tf = dynamic_cast<TotalField*>(source.get());
        if (tf == nullptr) continue;

        for (int i = 0; i < ndofs; ++i) {
            const double sign = tfsf_sign_[i];
            for (auto ft : { E, H }) {
                for (auto d : { X, Y, Z }) {
                    cached_tfsf_fields_[ft][d][i] +=
                        tf->eval(dof_coords_[i], time, static_cast<FieldType>(ft), static_cast<Direction>(d)) * sign;
                }
            }
        }
    }
}

}