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

        std::function<double(const GaussianInitialField::Position&)> f = 
            std::bind(&GaussianInitialField::eval1D, &source, std::placeholders::_1);

        switch (source.getFieldType()) {
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
        std::function<double(const GaussianInitialField::Position&)> f = 0;

        switch (fes_.GetMesh()->Dimension()) {
        case 2:
            f = std::bind(&GaussianInitialField::eval2D, &source, std::placeholders::_1);
            break;
        case 3:
            f = std::bind(&GaussianInitialField::eval3D, &source, std::placeholders::_1);
            break;
        default:
            throw std::exception("Incorrect Dimension for setFields3D");
        }

        switch (source.getFieldType()) {
        case FieldType::E:
            fields.E[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
            break;
        case FieldType::H:
            fields.H[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
            break;
        }
    }
}


}