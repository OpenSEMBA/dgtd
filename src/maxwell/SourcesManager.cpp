#include "SourcesManager.h"

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(const Sources& srcs, mfem::FiniteElementSpace& fes) :
    fes_{ fes }
{
    for (const auto& src : srcs) {
        sources.push_back(src->clone());
    }
}

void SourcesManager::setInitialFields(Fields& fields)
{
    for (const auto& source : sources) {
        if (dynamic_cast<InitialField*>(source.get())) {
            auto initialField{ dynamic_cast<InitialField*>(source.get()) };
            for (FieldType ft: {E, H}) {
                for (auto x : { X, Y, Z }) {
                    std::function<double(const Source::Position&, Source::Time)> f = 0;
                    f = std::bind(
                        &InitialField::eval, initialField, 
                        std::placeholders::_1, std::placeholders::_2, ft, x
                    );
                    FunctionCoefficient fc(f);
                    fields.get(ft, x).ProjectCoefficient(fc);
                }
            }
        }
    }
}

SourcesManager::TimeVarOperators SourcesManager::evalTimeVarField(const double time)
{
    SourcesManager::TimeVarOperators res;
    for (auto ft : { E, H }) {
        for (auto d : { X, Y, Z }) {
            res[ft][d].SetSpace(&fes_);
            res[ft][d] = 0.0;
        }
    }
    for (const auto& source : sources) {
        std::function<double(const Source::Position&, Source::Time)> f = 0;
        if (dynamic_cast<Planewave*>(source.get())) {
            auto timeVarField{ dynamic_cast<Planewave*>(source.get()) };
            f = std::bind(
               &Planewave::eval,
                timeVarField,
                std::placeholders::_1,
                std::placeholders::_2
            );
            FunctionCoefficient func(f);
            func.SetTime(time);
            for (auto d : { X, Y, Z }) {
                res[timeVarField->fieldType_][d].ProjectCoefficient(func);
                res[timeVarField->fieldType_][d] *= timeVarField->polarization_[d];
                res[timeVarField->fieldType_][d] *= 0.5;
            }
        }
    }
    return res;
}

}