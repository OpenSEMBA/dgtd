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
// File: simulation.h
#ifndef SOLVERDGTD_H_
#define SOLVERDGTD_H_

#include "../../exporter/gid/Exporter.h"
#include "exporter/Output.h"
#include "Options.h"
#include "parser/gid/Parser.h"
#include "physicalModel/Group.h"
#ifdef USE_MPI
    #include "CommMPI.h"
    #include "../output/OutputCommGiD.h"
    #include "../output/OutputComm.h"
#else
    #include "../dgtd/core/CommNone.h"
#endif
#include "../dgtd/integrator/IntegratorLSERK.h"
#include "../dgtd/integrator/IntegratorLF2.h"
#include "../dgtd/integrator/IntegratorLF2Full.h"
#include "../dgtd/integrator/IntegratorVerlet.h"
#include "solver/Solver.h"
#include "dg/DGExplicit.h"

class Solver {
public:
    Solver(SmbData*);
    ~Solver();
    bool run();
    bool canRun() const;
private:
    Comm *comm_;
    Integrator *integrator_;
    DG *dg_;
    Exporter* exporter_;
    const OptionsSolverDGTD* options_;
    Integrator* initIntegrator(
            const Mesh::Volume* mesh,
            const PMGroup* pMGroup,
            const OptionsSolverDGTD* args);
    Comm* initMPI();
};

#endif
