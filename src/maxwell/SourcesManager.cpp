#include "SourcesManager.h"

namespace maxwell {

using namespace mfem;

SourcesManager::SourcesManager(Sources srcs, const mfem::FiniteElementSpace& fes) :
	sources{srcs},
	fes_{fes}
{}

void SourcesManager::setFields1D(Fields& fields)
{
    for (const auto& source : sources) {

        switch (source.get()->initialFT) {
        case InitialFieldType::Gaussian:
            std::function<double(const Source::Position&)> f = std::bind(&Source::eval1D, &source, std::placeholders::_1);
            break;
        case InitialFieldType::PlanarSinusoidal:
            std::function<double(const Source::Position&)> f = 0;
            f = std::bind(&PlanarSinusoidalInitialField::eval1D, &source, std::placeholders::_1);
            break;
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

void SourcesManager::setGaussianSource(std::unique_ptr<Source> source) 
{
    std::function<double(const GaussianInitialField::Position&)> f = std::bind(&GaussianInitialField::eval1D, &GaussianInitialField, std::placeholders::_1);
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