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
    //for (const auto& source : sources) {
    //    std::function<double(const GaussianInitialField::Position&)> f = 0;

    //    switch (fes_.GetMesh()->Dimension()) {
    //    case 2:
    //        f = std::bind(&GaussianInitialField::eval2D, &source, std::placeholders::_1);
    //        break;
    //    case 3:
    //        f = std::bind(&GaussianInitialField::eval2D, &source, std::placeholders::_1);
    //        break;
    //    default:
    //        throw std::exception("Incorrect Dimension for setFields3D");
    //    }

    //    switch (source.getFieldType()) {
    //    case FieldType::E:
    //        fields.E[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
    //        break;
    //    case FieldType::H:
    //        fields.H[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
    //        break;
    //    }
    //}
}


}