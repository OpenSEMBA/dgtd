#include "SourcesManager.h"

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

}