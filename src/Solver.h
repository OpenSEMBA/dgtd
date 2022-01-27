#pragma once

#include "mfem.hpp"

namespace Maxwell {

class Solver {
public:
    struct Options {
        int order = 3;
        double t_final = 10.0;
        double dt = 0.01;
        int vis_steps = 5;
        int precision = 8;
    };

	Solver(const Options&, const mfem::Mesh&);
    
    void setInitialField(std::function<double(const mfem::Vector&)>);
    
    mfem::Mesh& getMesh() { return mesh_;  }

    void run();

private:
    
    Options opts_;

    std::unique_ptr<mfem::FiniteElementSpace> fes_;

    mfem::Mesh mesh_;

    std::unique_ptr<mfem::BilinearForm> MInv_;
    std::unique_ptr<mfem::BilinearForm> Kx_;
    std::unique_ptr<mfem::BilinearForm> Ky_;

    mfem::GridFunction ez_, hx_, hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;
    
    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    //std::unique_ptr<mfem::FiniteElementSpace> buildFiniteElementSpace() const;
    void initializeBilinearForms();
    void buildDomainAndFaceIntegrators();
    void buildBilinearForms();
    void Solver::collectParaviewData();
};

}

