#pragma once

#include "Types.h"
#include "Fields.h"
#include "ProbesManager.h"
#include "SourcesManager.h"
#include "SolverOptions.h"
#include "evolution/Evolution3D.h"
#include "evolution/Evolution3D_Spectral.h"

namespace maxwell {

struct ProblemDescription {
    Model model;
    Probes probes;
    Sources sources;
};

class Solver {
public:
    using Vector = mfem::Vector;
    using Position = Vector;
    using GridFunction = mfem::GridFunction;
    using ODESolver = mfem::ODESolver;
    
    Solver(const ProblemDescription&, const SolverOptions& = SolverOptions());
    Solver(const Model&, const Probes&, const Sources&, const SolverOptions& = SolverOptions());
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    const Fields& getFields() const { return fields_; };
    const PointProbe& getPointProbe(const std::size_t probe) const;

    double getTime() const { return time_; }
    double getTimeStep();

    const TimeDependentOperator* getFEEvol() const { return maxwellEvol_.get(); }

    void run();
    void step();
private:
    SolverOptions opts_;
    Model model_;
    mfem::DG_FECollection fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;
    Fields fields_;
    
    SourcesManager sourcesManager_;
    ProbesManager probesManager_;
    
    double time_;
    double dt_;
    std::unique_ptr<ODESolver> odeSolver_{ std::make_unique<mfem::RK4Solver>() };
    
    std::unique_ptr<mfem::TimeDependentOperator> maxwellEvol_;

    void checkOptionsAreValid(const SolverOptions&) const;
    void initializeFieldsFromSources();
};
}