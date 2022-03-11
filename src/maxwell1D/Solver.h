#pragma once

#include "mfem.hpp"

#include "FE_Evolution.h"

namespace Maxwell1D {

class Solver {
public:
    enum class TimeIntegrator {
        Leapfrog,
        RK4
    };

    typedef double ElectricField;
    typedef mfem::Vector Position;
    typedef std::size_t Direction;
    typedef std::size_t FieldType;

    const Direction X = 0;

    const FieldType Electric = 0;
    const FieldType Magnetic = 1;

    struct Options {
        int order = 2;
        double dt = 1e-4;
        double t_final = 1000*dt;
        int vis_steps = 100;
        int precision = 8;
        bool paraview = false;
        bool glvis = false;
        //TimeIntegrator timeIntegrator = Leapfrog;
    };

    Solver(const Options&, const mfem::Mesh&);

    void setInitialElectricField(std::function<ElectricField(const Position&)>);

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

};

}