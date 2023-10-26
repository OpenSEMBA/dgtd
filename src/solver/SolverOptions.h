#pragma once

#include "evolution/EvolutionOptions.h"

namespace maxwell {

struct SolverOptions {
    double timeStep = 0.0;
    double finalTime = 2.0;
    double cfl = 0.8;

    EvolutionOptions evolution;
    
    SolverOptions& setTimeStep(double t) 
    {
        timeStep = t;
        return *this;
    };

    SolverOptions& setFinalTime(double t) 
    { 
        finalTime = t; 
        return *this; 
    };

    SolverOptions& setCentered() 
    {
        evolution.fluxType = FluxType::Centered;
        return *this;
    };

    SolverOptions& setCFL(double cfl) 
    {
        cfl = cfl;
        return *this;
    }

    SolverOptions& setOrder(int orderIn) 
    {
        evolution.order = orderIn;
        return *this;
    }

    SolverOptions& setSpectralEO(bool spectral = true) {
        evolution.spectral = spectral;
        return *this;
    }
    //SolverOptions& setBasis(BasisType bst) {
    //    basis = bst;
    //    return *this;
    //}
};

}