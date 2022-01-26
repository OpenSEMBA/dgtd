#pragma once

#include "mfem.hpp"

namespace Maxwell {

class Solver {
public:

    struct Options {
        int order = 3;
        const char* device_config = "cpu";
        double t_final = 10.0;
        double dt = 0.01;

        bool paraview = false;
        int vis_steps = 5;

        int precision = 8;
    };

	Solver(const Options&, mfem::Mesh&);
    
    void run();

protected:
    /// Initial condition. 
    static double u0_function(const mfem::Vector& x)
    {

        // map to the reference [-1,1] domain
        mfem::Vector X(2);
        for (size_t i = 0; i < 2; i++) {
            double center = (meshBoundingBoxMin[i] + meshBoundingBoxMax[i]) * 0.5;
            X[i] = 2 * (x[i] - center) / (meshBoundingBoxMax[0] - meshBoundingBoxMin[0]);
        }

        return exp(-10. * (pow(X[0], 2) + pow(X[1], 2)));
        //return exp(-10. * pow(X[0], 2));
    }

private:
    static mfem::Vector meshBoundingBoxMin, meshBoundingBoxMax;

    Options opts_;

    std::unique_ptr<mfem::FiniteElementSpace> fes_;

    std::unique_ptr<mfem::BilinearForm> MInv_;
    std::unique_ptr<mfem::BilinearForm> Kx_;
    std::unique_ptr<mfem::BilinearForm> Ky_;

    mfem::GridFunction ez_, hx_, hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;
    
    void initializeFiniteElementSpace(mfem::Mesh& mesh);
    void initializeBilinearForms();
    void buildDomainAndFaceIntegrators();
    void buildBilinearForms();
};

}

