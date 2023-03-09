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
        std::function<double(const Source::Position&, Source::Time)> f = 0;
        if (dynamic_cast<InitialField*>(source.get())) {
            auto initialField{ dynamic_cast<InitialField*>(source.get()) };
            f = std::bind(
                &InitialField::eval,
                initialField,
                std::placeholders::_1,
                std::placeholders::_2
            );
            switch (fes_.GetMesh()->Dimension()) {
            case 1:
                switch (initialField->fieldType) {
                case FieldType::E:
                    fields.E1D.ProjectCoefficient(FunctionCoefficient(f));
                    break;
                case FieldType::H:
                    fields.H1D.ProjectCoefficient(FunctionCoefficient(f));
                    break;
                }
                break;
            default:
                switch (initialField->fieldType) {
                case FieldType::E:
                    for (auto x : { X, Y, Z }) {
                        fields.E[x].ProjectCoefficient(FunctionCoefficient(f));
                        fields.E[x] *= initialField->polarization[x];
                    }
                    break;
                case FieldType::H:
                    for (auto x : { X, Y, Z }) {
                        fields.H[x].ProjectCoefficient(FunctionCoefficient(f));
                        fields.H[x] *= initialField->polarization[x];
                    }
                    break;
                }
                break;
            }
        }
    }
}

SourcesManager::TimeVarOperators SourcesManager::evalTimeVarField(const double time)
{
    SourcesManager::TimeVarOperators res;
    for (auto d : { X, Y, Z }) {
        res[d].SetSpace(&fes_);
    }
    for (const auto& source : sources) {
        std::function<double(const Source::Position&, Source::Time)> f = 0;
        if (dynamic_cast<TimeVaryingField*>(source.get())) {
            auto timeVarField{ dynamic_cast<TimeVaryingField*>(source.get()) };
            f = std::bind(
               &TimeVaryingField::eval,
                timeVarField,
                std::placeholders::_1,
                std::placeholders::_2
            );
            FunctionCoefficient func(f);
            func.SetTime(time);
            for (auto f : { E, H }) {
                for (auto d : { X, Y, Z }) {
                    res[d].ProjectCoefficient(func);
                    res[d] *= timeVarField->polarization[d];
                }
            }
        }
    }
    return res;
}

}