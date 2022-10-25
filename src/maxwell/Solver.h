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

    //const MaxwellEvolution1D& getFEEvol1D() const { return maxwellEvol1D_; }
    const MaxwellEvolution2D& getFEEvol2D() const { return maxwellEvol2D_; }
    const MaxwellEvolution3D& getFEEvol3D() const { return maxwellEvol3D_; }

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
    //MaxwellEvolution1D maxwellEvol1D_;
    MaxwellEvolution2D maxwellEvol2D_;
    MaxwellEvolution3D maxwellEvol3D_;

    void checkOptionsAreValid(const SolverOptions&);

    const double Solver::calculateTimeStep() const;

    void Solver::initializeFieldsFromSources();
};
}