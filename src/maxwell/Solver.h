#pragma once

#include "Types.h"
#include "Fields.h"
#include "ProbesManager.h"
#include "SourcesManager.h"
#include "SolverOptions.h"
#include "MaxwellEvolution3D.h"
#include "MaxwellEvolution2D.h"
#include "MaxwellEvolution1D.h"

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
    const PointsProbe& getPointsProbe(const std::size_t probe) const;

    const TimeDependentOperator* getFEEvol() const { return maxwellEvol_.get(); }

    void run();

private:
    SolverOptions opts_;
    Model model_;
    mfem::DG_FECollection fec_;
    mfem::FiniteElementSpace fes_;
    Fields fields_;
    
    SourcesManager sourcesManager_;
    ProbesManager probesManager_;
    
    double time_;
    std::unique_ptr<ODESolver> odeSolver_{ std::make_unique<mfem::RK4Solver>() };
    
    std::unique_ptr<mfem::TimeDependentOperator> maxwellEvol_;

    void checkOptionsAreValid(const SolverOptions&);

    const double Solver::calculateTimeStep() const;

    void Solver::initializeFieldsFromSources();
};
}