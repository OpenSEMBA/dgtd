#pragma once

#include "Types.h"
#include "evolution/EvolutionOptions.h"

namespace maxwell {

struct SolverOptions {
    int order = 2;
    double dt = 0.0;
    double t_final = 2.0;
    double CFL = 0.8;

    EvolutionOptions evolutionOperatorOptions;
    
    SolverOptions& setTimeStep(double t) 
    {
        dt = t;
        return *this;
    };

    SolverOptions& setFinalTime(double t) 
    { 
        t_final = t; 
        return *this; 
    };

    SolverOptions& setCentered() 
    {
        evolutionOperatorOptions.fluxType = FluxType::Centered;
        return *this;
    };

    SolverOptions& setCFL(double cfl) 
    {
        CFL = cfl;
        return *this;
    }

    SolverOptions& setOrder(int orderIn) 
    {
        order = orderIn;
        return *this;
    }

    SolverOptions& setSpectralEO(bool spectral = true) {
        evolutionOperatorOptions.spectral = spectral;
        return *this;
    }
    //SolverOptions& setBasis(BasisType bst) {
    //    basis = bst;
    //    return *this;
    //}
};

}