#pragma once

#include "evolution/EvolutionOptions.h"

namespace maxwell {

// SolverOptions.h (or wherever ODEType lives)
enum ODEType : size_t {
    RK4              = 0, // [explicit] 4th order

    BackwardEuler    = 1, // [implicit] 1st order (dissipative)
    Trapezoidal      = 2, // [implicit] Crank–Nicolson, A-stable, low dissipation
    ImplicitMidpoint = 3, // [implicit] symplectic/energy-conserving for linear problems
    SDIRK33          = 4, // [implicit] A-stable, 3rd order
    SDIRK23          = 5, // [implicit] L-stable (use gamma_opt=2), ~2nd/3rd order
    SDIRK34          = 6  // [implicit] SDIRK(3/4) family, A-stable
};


struct SolverOptions {
    double timeStep = 0.0;
    double finalTime = 2.0;
    double cfl = 1.0;
    int basisType = mfem::BasisType::GaussLobatto;
    size_t odeType = ODEType::RK4;

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

    SolverOptions& setODEType(size_t ode_type) {
        odeType = ode_type;
        return *this;
    }

};

}