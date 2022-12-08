#pragma once

#include "ProblemDescription.h"

//#include "dg/Explicit.h"
//#include "communications/None.h"

namespace SEMBA::dgtd {


class Cudg3d {
public:
    struct Options {
        enum class TimeIntegrator {
            lserk4, verlet, lf2, lf2full
        };

        Math::Real upwinding{ 1.0 };
        TimeIntegrator timeIntegrator{ TimeIntegrator::lserk4 };
        bool useLTS{ true };
        std::size_t growSmallerTiers = 0;
        std::size_t maxNumberOfTiers = 0;
        bool useMaxStageSizeForLTS = false;
        bool PMLConstantConductivityProfile = false;
        Math::Real PMLConductivity = 0.0;
        Math::Real finalTime = 30e-9;
    };

    Cudg3d(const UnstructuredProblemDescription&, const Options&);
    void run();

private:
    Options options_;
//    Communications::Comm *comm_;
//    Integrator *integrator_;
//    DG *dg_;

//    Integrator* initIntegrator(
//            const Mesh::Volume* mesh,
//            const PMGroup* pMGroup,
//            const Options* args);
//    Communications::Comm* initMPI();
};

}