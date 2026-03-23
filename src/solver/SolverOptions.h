#pragma once

#include "evolution/EvolutionOptions.h"
#include "components/Material.h"

namespace maxwell {

enum ode_type : size_t {
    RK4              = 0, // [explicit] 4th order

    BackwardEuler    = 1, // [implicit] 1st order (dissipative)
    Trapezoidal      = 2, // [implicit] Crank–Nicolson, A-stable, low dissipation
    ImplicitMidpoint = 3, // [implicit] symplectic/energy-conserving for linear problems
    SDIRK33          = 4, // [implicit] A-stable, 3rd order
    SDIRK23          = 5, // [implicit] L-stable (use gamma_opt=2), ~2nd/3rd order
    SDIRK34          = 6  // [implicit] SDIRK(3/4) family, A-stable
};

struct SolverOptions {
    double time_step = 0.0;
    double final_time = 2.0;
    double cfl = 1.0;
    int basis_type = mfem::BasisType::GaussLobatto;
    size_t ode_type = ode_type::RK4;
    bool is_sgbc_solver = false;  // If true, skip statistics writing (SGBC sub-solver)

    EvolutionOptions evolution;
    
    SolverOptions& setTimeStep(double t) 
    {
        time_step = t;
        return *this;
    };

    SolverOptions& setFinalTime(double t) 
    { 
        final_time = t; 
        return *this; 
    };

    SolverOptions& setUpwindAlpha(double alpha) 
    {
        evolution.alpha = alpha;
        return *this;
    };

    SolverOptions& setCFL(double c) 
    {
        cfl = c;
        return *this;
    }

    SolverOptions& setOrder(int o) 
    {
        evolution.order = o;
        return *this;
    }

    SolverOptions& setSpectralEO(bool spectral = true) {
        evolution.spectral = spectral;
        return *this;
    }

    SolverOptions& setExportEO(bool export_op = true) {
        evolution.export_evolution_operator = export_op;
        return *this;
    }

    SolverOptions& setEvolutionOperator(EvolutionOperatorType oper) {
        evolution.op = oper;
        return *this;
    }    

    SolverOptions& setBasisType(int bt = mfem::BasisType::GaussLobatto) {
        basis_type = bt;
        return *this;
    }

    SolverOptions& setTFSFCutoffTime(double tfsf_ft) {
        evolution.tfsf_cutoff_time = tfsf_ft;
        return *this;
    }

    SolverOptions& setODEType(size_t type) {
        ode_type = type;
        return *this;
    }

    SolverOptions& setIsSGBCSolver(bool is_sgbc = true) {
        is_sgbc_solver = is_sgbc;
        evolution.is_sgbc_solver = is_sgbc;
        return *this;
    }

};

}