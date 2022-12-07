#pragma once


#include "core/physicalModel/Group.h"
#include "core/solver/Solver.h"

#include "Options.h"
#include "dg/Explicit.h"
//#include "communications/None.h"
//#include "mesh/Volume.h"

namespace SEMBA::dgtd {

class Solver {
public:
    Solver(Data* raw);
    virtual ~Solver() = default;

    void run();
private:
//    Communications::Comm *comm_;
//    Integrator *integrator_;
//    DG *dg_;
//    Exporter* exporter_;
    const Options* options_;

//    Integrator* initIntegrator(
//            const Mesh::Volume* mesh,
//            const PMGroup* pMGroup,
//            const Options* args);
//    Communications::Comm* initMPI();
};

}