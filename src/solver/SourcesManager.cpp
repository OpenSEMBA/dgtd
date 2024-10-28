#include "SourcesManager.h"

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(const Sources& srcs, mfem::FiniteElementSpace& fes, Fields& fields) :
    sources{ srcs }, 
    fes_{ fes }
{
    setInitialFields(fields);
}

void SourcesManager::setInitialFields(Fields& fields)
{
    for (const auto& source : sources) {
        auto src{ dynamic_cast<InitialField*>(source.get()) };
        if (src == nullptr) {
            continue;
        }
        for (FieldType ft: {E, H}) {
            for (auto x : { X, Y, Z }) {
                std::function<double(const Source::Position&, Source::Time)> f = 0;
                f = std::bind(
                    &InitialField::eval, src, 
                    std::placeholders::_1, std::placeholders::_2, ft, x
                );
                FunctionCoefficient fc(f);
                GridFunction gf(fields.get(ft, x).FESpace());
                gf.ProjectCoefficient(fc);
                fields.get(ft, x) += gf;
            }
        }
    }
}

FieldGridFuncs SourcesManager::evalTimeVarField(const Time time)
{
    std::array<std::array<GridFunction, 3>, 2> res;
    for (const auto& source : sources) {
        auto pw = dynamic_cast<Planewave*>(source.get());
        if (pw == nullptr) {
            continue;
        }
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                std::function<double(const Source::Position&, Source::Time)> f = 0;
                f = std::bind(&Planewave::eval, pw, 
                    std::placeholders::_1, std::placeholders::_2, ft, d);
                FunctionCoefficient func(f);
                func.SetTime(time);
                
                res[ft][d].SetSpace(&fes_);
                res[ft][d].ProjectCoefficient(func);
            }
        }
    }
    return res;
}

FieldGridFuncs SourcesManager::evalTimeVarField(const Time time, FiniteElementSpace* fes)
{
    std::array<std::array<GridFunction, 3>, 2> res;
    for (const auto& source : sources) {
        auto pw = dynamic_cast<Planewave*>(source.get());
        if (pw == nullptr) {
            continue;
        }
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                std::function<double(const Source::Position&, Source::Time)> f = 0;
                f = std::bind(&Planewave::eval, pw,
                    std::placeholders::_1, std::placeholders::_2, ft, d);
                FunctionCoefficient func(f);
                func.SetTime(time);
                res[ft][d].SetSpace(fes);
                res[ft][d].ProjectCoefficient(func);
            }
        }
    }
    return res;
}

FieldGridFuncs SourcesManager::evalTimeVarField(const Time time, bool is_tf)
{
    std::array<std::array<GridFunction, 3>, 2> res;
    for (const auto& source : sources) {
        auto pw = dynamic_cast<Planewave*>(source.get());
        if (pw == nullptr) {
            continue;
        }
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                std::function<double(const Source::Position&, Source::Time)> f = 0;
                f = std::bind(&Planewave::eval, pw,
                    std::placeholders::_1, std::placeholders::_2, ft, d);
                FunctionCoefficient func(f);
                func.SetTime(time);

                switch (is_tf) {
                case true:
                    res[ft][d].SetSpace(tf_fes_.get());
                    break;
                case false:
                    res[ft][d].SetSpace(sf_fes_.get());
                    break;
                }
                res[ft][d].ProjectCoefficient(func);
            }
        }
    }
    return res;
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
                    gfs[f][d][dofs[i]] = 0.0;
                }
            }
        }
    }
}

void SourcesManager::initTFSFPreReqs(const Mesh& m, const Array<int>& marker)
{
    initTFSFSubMesher(m, marker);
    initTFSFSpaces();
}

void SourcesManager::initTFSFSubMesher(const Mesh& m, const Array<int>& marker)
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