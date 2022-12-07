#include "Solver.h"

namespace SEMBA {
namespace Cudg3d {
namespace Solver {

Solver::Solver(Data* raw) {
    // Smb data adaptation and validation.
    //    Mesh::Volume mesh(*raw->mesh->castTo<MeshUnstructured>());
    Data smb;
    //    AdapterDGTD(*raw).convert(smb); TODO Adapt OutRqs on Volumes.
    // Time integrator initialization.
    //    options_ = smb.solverOptions->castTo<OptionsSolverDGTD>();
    //    integrator_ = initIntegrator(&mesh, smb.pMGroup, options_);
    //    integrator_->partitionate(&mesh, comm_);
    // Spatial discretization.
    //    dg_ = new DGExplicit(mesh, *smb.physicalModels, *smb.sources, *options_, comm_);
    //    dg_->setOutputRequests(smb.outputRequests);
    //    integrator_->setSolver(dg_);
    // Exporter initialization.
    //    cout << " - Initializing exporter... " << flush;
    //    const string outputFilename = smb.getOutputFilename();
    //    exporter_ = new ExporterGiD(smb, outputFilename, outputs_);
    //    cout << "[OK]" << endl;

}

bool Solver::run() {
    Math::Real time = 0.0;
    while (time < options_->getFinalTime()) {
        //        dg_->update(outputs_);
        //        exporter_->process(time, outputs_);
    }
    return true;
}

//Integrator* Solver::initIntegrator(
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
}
