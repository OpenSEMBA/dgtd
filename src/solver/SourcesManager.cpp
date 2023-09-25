#include "SourcesManager.h"

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(const Sources& srcs, mfem::FiniteElementSpace& fes) :
    sources{ srcs }, 
    fes_{ fes }, 
    tf_fes_{ FiniteElementSpace{&tf_mesh_,fes_.FEColl()} }, 
    sf_fes_{ FiniteElementSpace{&sf_mesh_,fes_.FEColl()} }
{
    tf_mesh_

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

void SourcesManager::initTFSFmeshes(const std::pair<SubMesh, SubMesh>& ms)
{
    auto tf_mesh{ ms.first };
    auto sf_mesh{ ms.second };
    tf_mesh_ = std::move(tf_mesh);
    sf_mesh_ = std::move(sf_mesh);

}

std::array<std::array<GridFunction, 3>, 2> SourcesManager::evalTimeVarField(const Time time)
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

std::array<std::array<GridFunction, 3>, 2> SourcesManager::evalTimeVarField(const Time time, bool is_tf)
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
                    res[ft][d].SetSpace(&tf_fes_);
                    break;
                case false:
                    res[ft][d].SetSpace(&sf_fes_);
                    break;
                }
                res[ft][d].ProjectCoefficient(func);
            }
        }
    }
    return res;
}

}