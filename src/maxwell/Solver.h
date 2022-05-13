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
    using EMFieldByVDIM = std::array<std::array<double, 3>, 3>;
    
    struct Options {
        int order = 2;
        double dt = 1e-3;
        double t_final = 1.0;
        FiniteElementEvolutionNoCond::Options evolutionOperatorOptions;
    };

    Solver(const Model&, const Probes&, const Sources&, const Options&);

    void setInitialField();
    const GridFunction& getFieldInDirection(const FieldType&, const Direction&) const;
    const Vector& getMaterialProperties(const Material&) const;

    mfem::Mesh& getMesh() { return mesh_; }
    std::vector<std::pair<double, EMFieldByVDIM>>& getFieldAtPoint() { return timeField_; }

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

    std::unique_ptr<FiniteElementEvolutionNoCond> maxwellEvol_;

    Vector sol_;

    std::array<GridFunction, 3> E_, H_;

    Array<int> elemIds_;
    IntegrationPointsSet integPointSet_;
    FieldType fieldToExtract_;
    double timeRecord_;
    EMFieldByVDIM fieldRecord_;
    std::vector<std::pair<double, EMFieldByVDIM>> timeField_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    socketstream sout_;

    void checkOptionsAreValid(const Options&);

    std::pair<Array<int>,Array<IntegrationPoint>>& Solver::buildElemAndIntegrationPointArrays(DenseMatrix& physPoints);
    const IntegrationPointsSet& Solver::buildIntegrationPointsSet(const Array<IntegrationPoint>& ipArray) const;
    const EMFieldByVDIM& saveFieldAtPoints(const FieldType&);

    void initializeParaviewData();
    //void initializeGLVISData();
    void storeInitialVisualizationValues();

};
}