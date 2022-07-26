#pragma once

#include "mfem.hpp"
#include "Material.h"
#include "FiniteElementEvolution.h"
#include "Types.h"
#include "Model.h"
#include "Probes.h"
#include "Sources.h"
#include "Options.h"
#include <array>

namespace maxwell {

class Solver {
public:

    using IntegrationPointsSet = std::vector<std::vector<IntegrationPoint>>;
    
    struct Options {
        int order = 2;
        double dt = 1e-3;
        double t_final = 1.0;
        FiniteElementEvolution::Options evolutionOperatorOptions;
    };

    Solver(const Model&, Probes&, const Sources&, const Options&);

    const GridFunction& getFieldInDirection(const FieldType&, const Direction&) const;
    const PointsProbe& getPointsProbe(const std::size_t probe) { return probes_.getPointsProbes().at(probe); }

    mfem::Mesh& getMesh() { return mesh_; }
    const mfem::Mesh& getConstMesh() const { return mesh_; }
    const std::unique_ptr<FiniteElementEvolution>& getFEEvol() const { return maxwellEvol_; }

    void run();

private:

    Model model_;
    Probes probes_;
    Sources sources_;
    Options opts_;
    
    mfem::Mesh& mesh_;

    std::unique_ptr<mfem::DG_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;

    std::unique_ptr<ODESolver> odeSolver_;

    mfem::Array<int> boundaryTDoF_;

    std::unique_ptr<FiniteElementEvolution> maxwellEvol_;

    Vector sol_;

    std::array<GridFunction, 3> E_, H_;

    std::vector<Array<int>> elemIds_;
    std::vector<IntegrationPointsSet> integPointSet_;
    double timeRecord_;
    std::vector<FieldFrame> fieldRecord_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    socketstream sout_;

    void checkOptionsAreValid(const Options&);

    void Solver::initializeSources();

    const std::pair<Array<int>, Array<IntegrationPoint>> buildElemAndIntegrationPointArrays(DenseMatrix& physPoints) const;
    const IntegrationPointsSet buildIntegrationPointsSet(const Array<IntegrationPoint>& ipArray) const;
    const std::vector<FieldFrame> saveFieldAtPointsForAllProbes();

    void initializeParaviewData();
    //void initializeGLVISData();
    void storeInitialVisualizationValues();

};
}