// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
#include "Solver.h"

namespace SEMBA {
namespace Cudg3d {

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

Solver::~Solver() {
//    delete exporter_;
//    delete dg_;
//    delete integrator_;
//    delete comm_;
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

bool Solver::canRun() const {
    throw logic_error("Cudg3d can run is not implemented.");
}

}
}
