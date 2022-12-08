#include "Cudg3d.h"

namespace SEMBA {
namespace dgtd {

Cudg3d::Cudg3d(const UnstructuredProblemDescription& raw, const Options& opts)
{
    // Smb data adaptation and validation.
    //    Mesh::Volume mesh(*raw->mesh->castTo<MeshUnstructured>());
    //    AdapterDGTD(*raw).convert(smb);
     
    dg_ = std::make_unique<dg::Evolution>(raw.model, raw.sources, opts.evolution);
    // Time integrator initialization.
    //    options_ = smb.solverOptions->castTo<OptionsSolverDGTD>();
    //    integrator_ = initIntegrator(&mesh, smb.pMGroup, options_);
    //    integrator_->partitionate(&mesh, comm_);
    //    integrator_->setSolver(dg_);
}

void Cudg3d::run() {
    Math::Real time = 0.0;
    while (time < options_.finalTime) {
        //        dg_->update(outputs_);
        //        exporter_->process(time, outputs_);
    }
}

//Integrator* Cudg3d::initIntegrator(
//        const Mesh::Volume* mesh,
//        const PMGroup* pMGroup,
//        const Options* arg) {
//    Integrator* res;
//    switch (arg->getTimeIntegrator()) {
//    case Options::TimeIntegrator::lserk4:
//        cout<< "- Initial   izing LSERK Integrator." << endl;
//        res = new IntegratorLSERK(*mesh, *pMGroup, arg);
//        break;
//    case Options::TimeIntegrator::lf2:
//        cout<< "- Initializing LF2 Integrator." << endl;
//        res = new IntegratorLF2(*mesh, *pMGroup, arg);
//        break;
//    case Options::TimeIntegrator::lf2full:
//        cout<< "- Initializing LF2Full Integrator." << endl;
//        res = new IntegratorLF2Full(*mesh, *pMGroup, arg);
//        break;
//    case Options::TimeIntegrator::verlet:
//        cout<< "- Initializing Verlet Integrator." << endl;
//        res = new IntegratorVerlet(*mesh, *pMGroup, arg);
//        break;
//    default:
//        throw logic_error("Undefined time integrator.");
//    }
//    return res;
//}

}
}
