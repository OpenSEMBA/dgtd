#include "SourcesManager.h"

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(const Sources& srcs, mfem::FiniteElementSpace& fes) :
    sources{srcs}, fes_{ fes }
{
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

std::array<std::array<GridFunction, 3>, 2> SourcesManager::evalTimeVarField(const double time)
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

}