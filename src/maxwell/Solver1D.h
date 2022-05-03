#pragma once

#include "mfem.hpp"
#include "Material.h"
#include "FiniteElementEvolution.h"
#include "Types.h"
#include "Model.h"
#include "Probes.h"
#include "Sources.h"
#include <array>

namespace maxwell {

class Solver1D {
public:
    using Position = mfem::Vector;
    using IntegrationPointsSet = std::vector<std::vector<IntegrationPoint>>;
    

    struct Options {
        int order = 2;
        double dt = 1e-3;
        double t_final = 1.0;
        FiniteElementEvolutionNoCond::Options evolutionOperatorOptions;
    };

    Solver1D(const Model&, const Probes&, const Sources&, const Options&);

    void setInitialField(const FieldType&, std::function<double(const Position&)>, const Direction&);
    const GridFunction& getFieldInDirection(const FieldType&, const Direction&) const;
    const Vector& getMaterialProperties(const Material&) const;

    mfem::Mesh& getMesh() { return mesh_; }
    Vector& getFieldAtPoint() { return timeField_; }

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

    FieldType fieldToExtract_;
    Vector timeRecord_;
    Vector fieldRecord_;
    Vector timeField_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    socketstream sout_;

    void checkOptionsAreValid(const Options&);

    std::pair<Array<int>,Array<IntegrationPoint>> Solver1D::buildIntegrationPointAndElemArrays(DenseMatrix& physPoints);
    const IntegrationPointsSet Solver1D::buildIntegrationPointsSet(const Array<IntegrationPoint>& ipArray) const;
    const std::array<std::array<double, 3>, 3> saveFieldAtPoints(DenseMatrix& physPoints, const FieldType&);

    void initializeParaviewData();
    //void initializeGLVISData();
    void storeInitialVisualizationValues();

};
}