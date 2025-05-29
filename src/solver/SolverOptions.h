#pragma once

#include "evolution/EvolutionOptions.h"

namespace maxwell {

struct SolverOptions {
    double timeStep = 0.0;
    double finalTime = 2.0;
    double cfl = 0.8;
    int basisType = mfem::BasisType::GaussLobatto;

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

    SolverOptions& setUpwindAlpha(double alpha) 
    {
        evolution.alpha = alpha;
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

    SolverOptions& setEvolutionOperator(EvolutionOperatorType oper) {
        evolution.op = oper;
        return *this;
    }    

    SolverOptions& setBasisType(int bt = mfem::BasisType::GaussLobatto) {
        basisType = bt;
        return *this;
    }
};

}