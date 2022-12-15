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

void SourcesManager::setFields1D(Fields& fields)
{
    for (const auto& source : sources) {

        std::function<double(const Source::Position&)> f = 0;
        if (dynamic_cast<GaussianInitialField*>(source.get())) {
            f = std::bind(
                &GaussianInitialField::eval1D, 
                dynamic_cast<GaussianInitialField*>(source.get()),
                std::placeholders::_1
            );
        }
        else if (dynamic_cast<PlanarSinusoidalInitialField*>(source.get())) {
            f = std::bind(
                &PlanarSinusoidalInitialField::eval1D, 
                dynamic_cast<PlanarSinusoidalInitialField*>(source.get()), 
                std::placeholders::_1
            );
        }
        else {
            throw std::runtime_error("Invalid source type.");
        }

        switch (source.get()->fieldType) {
        case FieldType::E:
            fields.E1D.ProjectCoefficient(FunctionCoefficient(f));
            break;
        case FieldType::H:
            fields.H1D.ProjectCoefficient(FunctionCoefficient(f));
            break;
        }
    }
}

void SourcesManager::setFields3D(Fields& fields)
{
    for (const auto& source : sources) {

        std::function<double(const Source::Position&)> f = 0;
        switch (fes_.GetMesh()->Dimension()) {
        case 2:
            if (dynamic_cast<GaussianInitialField*>(source.get())) {
                f = std::bind(
                    &GaussianInitialField::eval2D,
                    dynamic_cast<GaussianInitialField*>(source.get()),
                    std::placeholders::_1
                );
            }
            else if (dynamic_cast<PlanarSinusoidalInitialField*>(source.get())) {
                f = std::bind(
                    &PlanarSinusoidalInitialField::eval2D,
                    dynamic_cast<PlanarSinusoidalInitialField*>(source.get()),
                    std::placeholders::_1
                );
            }
            else {
                throw std::runtime_error("Invalid source type.");
            }
            break;
        case 3:
            if (dynamic_cast<GaussianInitialField*>(source.get())) {
                f = std::bind(
                    &GaussianInitialField::eval3D,
                    dynamic_cast<GaussianInitialField*>(source.get()),
                    std::placeholders::_1
                );
            }
            else if (dynamic_cast<PlanarSinusoidalInitialField*>(source.get())) {
                f = std::bind(
                    &PlanarSinusoidalInitialField::eval3D,
                    dynamic_cast<PlanarSinusoidalInitialField*>(source.get()),
                    std::placeholders::_1
                );
            }
            else {
                throw std::runtime_error("Invalid source type.");
            }
            break;
        default:
            throw std::exception("Incorrect Dimension for setFields3D");
        }

        switch (source.get()->fieldType) {
        case FieldType::E:
            fields.E[source.get()->direction].ProjectCoefficient(FunctionCoefficient(f));
            break;
        case FieldType::H:
            fields.H[source.get()->direction].ProjectCoefficient(FunctionCoefficient(f));
            break;
        }
    }
}


}