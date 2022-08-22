#pragma once

#include "Material.h"
#include "FiniteElementEvolution.h"
#include "ProbesManager.h"
#include "Types.h"
#include "SolverOptions.h"

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

    const GridFunction& getFieldInDirection(const FieldType&, const Direction&) const;
    const PointsProbe* getPointsProbe(const std::size_t probe) { return probesManager_.getPointsProbe(probe); }

    const FiniteElementEvolution& getFEEvol() const { return maxwellEvol_; }

    void run();

private:
    SolverOptions opts_;
    
    Model model_;
    Sources sources_;
    ProbesManager probesManager_;
    
    mfem::DG_FECollection fec_;
    mfem::FiniteElementSpace fes_;

    std::unique_ptr<ODESolver> odeSolver_{ std::make_unique<mfem::RK4Solver>() };

    FiniteElementEvolution maxwellEvol_;

    Vector sol_;
    std::array<GridFunction, 3> E_, H_;

    void checkOptionsAreValid(const SolverOptions&);

    void Solver::initializeFieldsFromSources();
};
}