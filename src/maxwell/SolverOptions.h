#pragma once

#include "FiniteElementEvolution.h"
#include "Types.h"

namespace maxwell {

struct SolverOptions {
    int order = 2;
    double dt = 1e-3;
    double t_final = 2.0;
    FiniteElementEvolution::Options evolutionOperatorOptions;

    SolverOptions& setFinalTime(double t) { 
        t_final = t; 
        return *this; 
    };
    SolverOptions& setCentered() {
        evolutionOperatorOptions.fluxType = FluxType::Centered;
        return *this;
    }
    SolverOptions& setStrongForm() {
        evolutionOperatorOptions.disForm = DisForm::Strong;
        return *this;
    }
};

}