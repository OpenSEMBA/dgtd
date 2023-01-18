#include "SourcesManager.h"

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(const Sources& srcs, const mfem::FiniteElementSpace& fes) :
	fes_{fes}
{
    for (const auto& src : srcs) {
        sources.push_back(src->clone());
    }
}

void SourcesManager::setFields(Fields& fields)
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
            default:
                switch (initialField->fieldType) {
                case FieldType::E:
                    fields.E[initialField->direction].ProjectCoefficient(FunctionCoefficient(f));
                    break;
                case FieldType::H:
                    fields.H[initialField->direction].ProjectCoefficient(FunctionCoefficient(f));
                    break;
                }
                break;
            }
        }
        else {
            throw std::runtime_error("Invalid source type.");
        }
        
    }
}

}