#pragma once

#include "Types.h"

namespace maxwell {

struct SolverOptions {
    int order = 2;
    double dt = 1e-3;
    double t_final = 2.0;
    double CFL = 0.9;
    MaxwellEvolOptions evolutionOperatorOptions;
    
    SolverOptions& setTimeStep(double t) {
        dt = t;
        return *this;
    };
    SolverOptions& setFinalTime(double t) { 
        t_final = t; 
        return *this; 
    };
    SolverOptions& setCentered() {
        evolutionOperatorOptions.fluxType = FluxType::Centered;
        return *this;
    };
    SolverOptions& setCFL(double cfl) {
        CFL = cfl;
        return *this;
    }
    SolverOptions& setOrder(int or) {
        order = or;
        return *this;
    }
};

}