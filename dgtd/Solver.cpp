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

Solver::Solver(SmbData* raw) {

    // Smb data adaptation and validation.
    Mesh::Volume mesh(*raw->mesh->castTo<MeshUnstructured>());
    SmbData smb;
//    AdapterDGTD(*raw).convert(smb); TODO Adapt OutRqs on Volumes.

    // Time integrator initialization.
    options_ = smb.solverOptions->castTo<OptionsSolverDGTD>();
    integrator_ = initIntegrator(&mesh, smb.pMGroup, options_);
    integrator_->partitionate(&mesh, comm_);

    // Spatial discretization.
    dg_ = new DGExplicit(mesh, *smb.pMGroup, *smb.emSources, *options_, comm_);
    dg_->setOutputRequests(smb.outputRequests);
    integrator_->setSolver(dg_);

    // Exporter initialization.
//    cout << " - Initializing exporter... " << flush;
//    const string outputFilename = smb.getOutputFilename();
//    exporter_ = new ExporterGiD(smb, outputFilename, outputs_);
//    cout << "[OK]" << endl;
}

Solver::~Solver() {
    delete exporter_;
    delete dg_;
    delete integrator_;
    delete comm_;
}

bool Solver::run() {
    Real tSum = 0.0;
    Real tRunning = 0.0;
    Real time = 0.0;
    const Real dt = integrator_->getMaxDT();
    assert(dt != 0.0);
    while (time < options_->getFinalTime()) {
//        dg_->update(outputs_);
//        exporter_->process(time, outputs_);
        Real initCPUTime = storeCPUTime();
        integrator_->timeIntegrate(time);
        tSum += storeCPUTime() - initCPUTime;
        time += dt;
        printTimeProfilingInfo(tSum, tRunning, time / dt,
                options_->getNumberOfTimeSteps());
    }
    return true;
}

Integrator* Solver::initIntegrator(
        const Mesh::Volume* mesh,
        const PMGroup* pMGroup,
        const OptionsSolverDGTD* arg) {
    Integrator* res;
    switch (arg->getTimeIntegrator()) {
    case OptionsSolverDGTD::lserk4:
        cout<< "- Initializing LSERK Integrator." << endl;
        res = new IntegratorLSERK(*mesh, *pMGroup, arg);
        break;
    case OptionsSolverDGTD::lf2:
        cout<< "- Initializing LF2 Integrator." << endl;
        res = new IntegratorLF2(*mesh, *pMGroup, arg);
        break;
    case OptionsSolverDGTD::lf2full:
        cout<< "- Initializing LF2Full Integrator." << endl;
        res = new IntegratorLF2Full(*mesh, *pMGroup, arg);
        break;
    case OptionsSolverDGTD::verlet:
        cout<< "- Initializing Verlet Integrator." << endl;
        res = new IntegratorVerlet(*mesh, *pMGroup, arg);
        break;
    default:
        throw Error("Undefined time integrator.");
    }
    return res;
}

bool Solver::canRun() const {
    throw ErrorNotImplemented("Can Run is not implemented.");
}
