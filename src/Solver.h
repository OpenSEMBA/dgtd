#pragma once

#include "mfem.hpp"

namespace Maxwell {

class Solver {
public:
    typedef double ElectricField;
    typedef mfem::Vector Position;

    typedef std::size_t Direction;

    const Direction X = 0;
    const Direction Y = 1;


    struct Options {
        int order = 3;
        double t_final = 10.0;
        double dt = 0.01;
        int vis_steps = 5;
        int precision = 8;
    };

	Solver(const Options&, const mfem::Mesh&);
    
    void setInitialElectricField(std::function<ElectricField(const Position&)>);
    
    mfem::Mesh& getMesh() { return mesh_;  }
    
    ElectricField getElectricFieldAtPosition(const Position&) const;

    void run();

private:
    
    Options opts_;

    std::unique_ptr<mfem::DG_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;
    std::unique_ptr<mfem::H1_FECollection> fecH1_;
    std::unique_ptr<mfem::FiniteElementSpace> fesH1_;

    mfem::Array<int> boundaryTDOFs_;

    mfem::Mesh mesh_;

    std::unique_ptr<mfem::BilinearForm> MInv_;
    std::unique_ptr<mfem::BilinearForm> Kx_;
    std::unique_ptr<mfem::BilinearForm> Ky_;

    mfem::GridFunction Ez_, hx_, Hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;
    
    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    mfem::Array<int> Solver::buildEssentialTrueDOF();
    std::unique_ptr<mfem::BilinearForm> buildMassMatrix() const;
    std::unique_ptr<mfem::BilinearForm> buildDerivativeOperator(const Direction&) const;

    void initializeParaviewData();
};

}

