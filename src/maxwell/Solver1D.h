#pragma once

#include "mfem.hpp"

#include "FE_Evolution.h"
#include "Types.h"

namespace maxwell {

class Solver1D {
public:
    typedef mfem::Vector Position;

    struct Options {
        int order = 2;
        double dt = 1e-3;
        double t_final = 1.0;
        int vis_steps = 1;
        int precision = 8;
        bool paraview = false;
        bool glvis = false;
        bool extractDataAtPoint = false;
        FieldType fieldToExtract = FieldType::Electric;
        IntegrationPoint integPoint;
        FE_Evolution::Options evolutionOperatorOptions;
    };

    Solver1D(const Options&, const mfem::Mesh&);

    void setInitialField(const FieldType&, std::function<double(const Position&)>);
    const GridFunction& getField(const FieldType&) const;

  
    mfem::Mesh& getMesh() { return mesh_; }
    Vector& getFieldAtPoint() { return timeField_; }

    void run();

private:

    Options opts_;
    mfem::Mesh mesh_;

    std::unique_ptr<mfem::DG_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;
    
    std::unique_ptr<ODESolver> odeSolver_;

    mfem::Array<int> boundaryTDoF_;

    std::unique_ptr<FE_Evolution> maxwellEvol_;
    
    Vector sol_;
    GridFunction E_, H_;

    IntegrationPoint integPoint_;
    FieldType fieldToExtract_;
    Vector timeRecord_;
    Vector fieldRecord_;
    Vector timeField_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    socketstream sout_;

    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    mfem::Array<int> buildEssentialTrueDOF();

    const IntegrationPoint setIntegrationPoint(const IntegrationPoint&) const;
    const int getElementIndexForPosition(const IntegrationPoint&) const;
    const Array<double> getVertexPositionInPhysicalCoords(const Array<int>& elementVertex) const;
    const IntegrationPoint getRelativePositionInElement(const int&, const IntegrationPoint&) const;
    const double saveFieldAtPoint(const IntegrationPoint&, const FieldType&) const;
    
    void initializeParaviewData();
    void initializeGLVISData();
    void storeInitialVisualizationValues();

};

}