#pragma once

#include "ProblemDescription.h"

#include "dg/Evolution.h"
#include "integrator/TimeIntegrator.h"

namespace SEMBA::dgtd {


class Cudg3d {
public:
    struct Options {        
        dg::Evolution::Options              evolution;
        integrator::TimeIntegrator::Options timeIntegrator;
    };

    Cudg3d(const UnstructuredProblemDescription&, const Options&);
    void run();

private:
    Options options_;
    std::unique_ptr<dg::Evolution> dg_;
    std::unique_ptr<integrator::TimeIntegrator> integrator_;
//    Communications::Comm *comm_;

    std::unique_ptr<integrator::TimeIntegrator> buildIntegrator(
            const dg::Evolution&,
            const integrator::TimeIntegrator::Options&);
//    Communications::Comm* initMPI();
};

}