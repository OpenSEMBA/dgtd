#pragma once

#include "mfem.hpp"

namespace Maxwell {

class Solver {
public:
    typedef double ElectricField;
    typedef mfem::Vector Position;

    typedef std::size_t Direction;
    typedef std::size_t FieldType;

    const Direction X = 0;
    const Direction Y = 1;

    const FieldType Electric = 0;
    const FieldType Magnetic = 1;


    struct Options {
        int order = 3;
        double t_final = 1000*dt;
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
    std::unique_ptr<mfem::BilinearForm> KxE_;
    std::unique_ptr<mfem::BilinearForm> KyE_;
    std::unique_ptr<mfem::BilinearForm> KxH_;
    std::unique_ptr<mfem::BilinearForm> KyH_;

    mfem::GridFunction Ez_, Hx_, Hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;
    
    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    mfem::Array<int> Solver::buildEssentialTrueDOF();
    std::unique_ptr<mfem::BilinearForm> buildInverseMassMatrix() const;
    std::unique_ptr<mfem::BilinearForm> buildDerivativeAndFluxOperator(const Direction&, const FieldType& ) const;

    void initializeParaviewData();
};

}

