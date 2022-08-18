#pragma once

#include "Material.h"
#include "FiniteElementEvolution.h"
#include "Probes.h"
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
    using IntegrationPoint = mfem::IntegrationPoint;
    using ODESolver = mfem::ODESolver;
    using IntegrationPointsSet = std::vector<std::vector<IntegrationPoint>>;
    
    Solver(const Model&, const Probes&, const Sources&, const SolverOptions& = SolverOptions());
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    const GridFunction& getFieldInDirection(const FieldType&, const Direction&) const;
    const PointsProbe& getPointsProbe(const std::size_t probe) { return probes_.getPointsProbes().at(probe); }

    const mfem::Mesh& getMesh() const { return mesh_; }
    const FiniteElementEvolution& getFEEvol() const { return maxwellEvol_; }

    void run();

private:
    SolverOptions opts_;
    Model model_;
    Sources sources_;

    // --- Probes Manager? --
    struct PointsProbeIntegrationPoints {
        mfem::Array<int> elemIds;
        IntegrationPointsSet integPointSet;
    };

    Probes probes_;
    std::map<const PointsProbe*, PointsProbeIntegrationPoints> probeIntegrationPoints_;

    PointsProbeIntegrationPoints buildElemAndIntegrationPointArrays(const PointsProbe&) const;
    const IntegrationPointsSet 
        buildIntegrationPointsSet(const mfem::Array<IntegrationPoint>& ipArray) const;
    FieldFrame getFieldForPointsProbe(const PointsProbe& p) const;
    // -----------------------


    mfem::Mesh& mesh_;

    mfem::DG_FECollection fec_;
    mfem::FiniteElementSpace fes_;

    std::unique_ptr<ODESolver> odeSolver_;

    mfem::Array<int> boundaryTDoF_;

    FiniteElementEvolution maxwellEvol_;

    Vector sol_;

    std::array<GridFunction, 3> E_, H_;
    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    mfem::socketstream sout_;

    void checkOptionsAreValid(const SolverOptions&);

    void Solver::initializeSources();


    void initializeParaviewData();
    void storeInitialVisualizationValues();

};
}