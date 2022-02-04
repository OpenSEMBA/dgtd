#pragma once

#include "mfem.hpp"

namespace Maxwell {

typedef std::size_t Direction;

const Direction X = 0;
const Direction Y = 1;


class Solver {
public:
    typedef double ElectricField;
    typedef mfem::Vector Position;

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

    mfem::Mesh mesh_;

    std::unique_ptr<mfem::LinearForm> inflowForm_;

    std::unique_ptr<mfem::BilinearForm> MInv_;
    std::unique_ptr<mfem::BilinearForm> Kx_;
    std::unique_ptr<mfem::BilinearForm> Ky_;

    mfem::GridFunction ez_, hx_, hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;
    
    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    std::unique_ptr<mfem::LinearForm> buildInflowForm() const;
    std::unique_ptr<mfem::BilinearForm> buildMassMatrix() const;
    std::unique_ptr<mfem::BilinearForm> buildDerivativeOperator(const Direction&) const;

    void initializeParaviewData();
};

}

