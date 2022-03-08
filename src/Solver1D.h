#pragma once

#include "mfem.hpp"

namespace Maxwell {


class Solver1D {
public:
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
    };

    Solver1D(const Options&, const mfem::Mesh&);

    void setInitialElectricField(std::function<ElectricField(const Position&)>);

    mfem::Mesh& getMesh() { return mesh_; }

    void run();
    void runODESolver();

private:

    Options opts_;

    std::unique_ptr<mfem::DG_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;
    std::unique_ptr<mfem::H1_FECollection> fecH1_;
    std::unique_ptr<mfem::FiniteElementSpace> fesH1_;

    mfem::Mesh mesh_;

    mfem::Array<int> boundaryTDoF_;

    std::unique_ptr<mfem::BilinearForm> MInv_;
    std::unique_ptr<mfem::BilinearForm> KxE_;
    std::unique_ptr<mfem::BilinearForm> KxH_;

    mfem::GridFunction Ez_, Hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;

    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    mfem::Array<int> Solver1D::buildEssentialTrueDOF();
    std::unique_ptr<mfem::BilinearForm> buildInverseMassMatrix() const;
    std::unique_ptr<mfem::BilinearForm> buildDerivativeAndFluxOperator(const Direction&, const FieldType&) const;
    void initializeParaviewData();

};

}