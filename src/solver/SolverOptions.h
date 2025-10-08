#pragma once

#include "evolution/EvolutionOptions.h"

namespace maxwell {

struct SolverOptions {
    double timeStep = 0.0;
    double finalTime = 2.0;
    double cfl = 1.0;
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

    SolverOptions& setExportEO(bool exportOP = true) {
        evolution.exportEvolutionOperator = exportOP;
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

    SolverOptions& setTFSFFinalTime(double tfsf_ft) {
        evolution.tfsfFinalTime = tfsf_ft;
        return *this;
    }

};

}