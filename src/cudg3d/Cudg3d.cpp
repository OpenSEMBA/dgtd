#include "Cudg3d.h"

#include "integrator/LSERK4.h"
#include "integrator/LF2.h"

namespace SEMBA {
namespace dgtd {

Cudg3d::Cudg3d(const UnstructuredProblemDescription& raw, const Options& opts)
{
    // Smb data adaptation and validation.
    //    Mesh::Volume mesh(*raw->mesh->castTo<MeshUnstructured>());
    //    AdapterDGTD(*raw).convert(smb);
    dg::VolumeModel adaptedModel{ raw.model };

    dg_ = std::make_unique<dg::Evolution>(adaptedModel, raw.sources, opts.evolution);   
    integrator_ = buildIntegrator(*dg_, opts.timeIntegrator);
}

void Cudg3d::run() {
    Math::Real time = 0.0;
    while (time < options_.timeIntegrator.finalTime) {
        //        dg_->update(outputs_);
        //        exporter_->process(time, outputs_);
    }
}

std::unique_ptr<integrator::TimeIntegrator> Cudg3d::buildIntegrator(
    const dg::Evolution& dg,
    const integrator::TimeIntegrator::Options& opts) 
{
    using timeIntegratorType = integrator::TimeIntegrator::Options::Type;
    switch (opts.timeIntegrator) {
    case timeIntegratorType::lserk4:
        return std::make_unique<integrator::LSERK4>(dg, opts);
    //case timeIntegratorType::lf2:
    //    return std::make_unique<integrator::LF2>(dg, opts);    
    default:
        throw std::logic_error("Undefined time integrator.");
    }
}

}
}
