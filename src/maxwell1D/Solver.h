#pragma once

#include "mfem.hpp"

#include "FE_Evolution.h"

namespace maxwell1D {

class Solver {
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
        FE_Evolution::Options evolutionOperatorOptions;
    };

    Solver(const Options&, const mfem::Mesh&);

    void setInitialField(const FieldType&, std::function<double(const Position&)>);
    const GridFunction& getField(const FieldType&) const;

    mfem::Mesh& getMesh() { return mesh_; }

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

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    socketstream sout_;

    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    mfem::Array<int> Solver::buildEssentialTrueDOF();
    
    void initializeParaviewData();
    void initializeGLVISData();
    void storeInitialVisualizationValues();

};

}