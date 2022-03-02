#pragma once

#include "mfem.hpp"

namespace Maxwell {

    typedef std::size_t Direction;

    const Direction X = 0;

    class Solver {
    public:
        typedef double ElectricField;
        typedef double Position;

        struct Options {
            int order = 2;
            double dt = 1e-4;
            double t_final = 1000*dt;
            int vis_steps = 100;
            int precision = 8;
        };

        Solver(const Options&, const mfem::Mesh&);

        void setInitialElectricField(std::function<ElectricField(const Position&)>);

        mfem::Mesh& getMesh() { return mesh_; }

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

        mfem::GridFunction ez_, hx_, hy_;

        std::unique_ptr<mfem::ParaViewDataCollection> pd_;

        void checkOptionsAreValid(const Options&, const mfem::Mesh&);
        mfem::Array<int> Solver::buildEssentialTrueDOF();
        std::unique_ptr<mfem::BilinearForm> buildMassMatrix() const;
        std::unique_ptr<mfem::BilinearForm> buildDerivativeOperator(const Direction&) const;

        void initializeParaviewData();
    };

}