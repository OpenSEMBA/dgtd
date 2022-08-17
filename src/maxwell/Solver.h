#pragma once

#include <mfem.hpp>
#include "Material.h"
#include "FiniteElementEvolution.h"
#include "Probes.h"
#include "Types.h"

namespace maxwell {

class Solver {
public:
    using Vector = mfem::Vector;
    using Position = Vector;
    using GridFunction = mfem::GridFunction;
    using IntegrationPoint = mfem::IntegrationPoint;
    using ODESolver = mfem::ODESolver;
    using IntegrationPointsSet = std::vector<std::vector<IntegrationPoint>>;
    
    struct Options {
        int order = 2;
        double dt = 1e-3;
        double t_final = 1.0;
        FiniteElementEvolution::Options evolutionOperatorOptions;
    };

    Solver(const Model&, Probes&, const Sources&, const Options&);
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    const GridFunction& getFieldInDirection(const FieldType&, const Direction&) const;
    const PointsProbe& getPointsProbe(const std::size_t probe) { return probes_.getPointsProbes().at(probe); }

    const mfem::Mesh& getMesh() const { return mesh_; }
    const FiniteElementEvolution& getFEEvol() const { return maxwellEvol_; }

    void run();

private:

    Model model_;
    Probes probes_;
    Sources sources_;
    Options opts_;
    
    mfem::Mesh& mesh_;

    mfem::DG_FECollection fec_;
    mfem::FiniteElementSpace fes_;

    std::unique_ptr<ODESolver> odeSolver_;

    mfem::Array<int> boundaryTDoF_;

    FiniteElementEvolution maxwellEvol_;

    Vector sol_;

    std::array<GridFunction, 3> E_, H_;

    std::vector<mfem::Array<int>> elemIds_;
    std::vector<IntegrationPointsSet> integPointSet_;
    double timeRecord_;
    std::vector<FieldFrame> fieldRecord_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    mfem::socketstream sout_;

    void checkOptionsAreValid(const Options&);

    void Solver::initializeSources();

    const std::pair<mfem::Array<int>, mfem::Array<IntegrationPoint>> 
        buildElemAndIntegrationPointArrays(mfem::DenseMatrix& physPoints) const;
    const IntegrationPointsSet 
        buildIntegrationPointsSet(const mfem::Array<IntegrationPoint>& ipArray) const;
    const std::vector<FieldFrame> saveFieldAtPointsForAllProbes();

    void initializeParaviewData();
    void storeInitialVisualizationValues();

};
}