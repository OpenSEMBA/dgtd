#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell {

Solver::Solver(const Options& opts, Mesh& mesh)
{
    opts_ = opts;

    Device device(opts_.device_config);
    mesh.GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax, std::max(opts_.order, 1));

    initializeFiniteElementSpace(mesh);

    initializeBilinearForms();

    buildDomainAndFaceIntegrators();

    assembleBilinearForms();

    finalizeBilinearForms();

    run();
}

void Solver::initializeFiniteElementSpace(Mesh& mesh)
{
    DG_FECollection fec(opts_.order, mesh.Dimension(), BasisType::GaussLobatto);
    fes_ = std::make_unique<FiniteElementSpace>(&mesh, &fec);
}

void Solver::initializeBilinearForms()
{
    MInv_ = std::make_unique<BilinearForm>(fes_.get());
    Kx_ = std::make_unique<BilinearForm>(fes_.get());
    Ky_ = std::make_unique<BilinearForm>(fes_.get());
}

void Solver::buildDomainAndFaceIntegrators()
{
    ConstantCoefficient zero(0.0), one(1.0), mOne(-1.0);
    Vector nxVec(2);  nxVec(0) = 1.0; nxVec(1) = 0.0;
    Vector nyVec(2);  nyVec(0) = 0.0; nyVec(1) = 1.0;
    Vector n1Vec(2);  n1Vec(0) = 1.0; n1Vec(1) = 1.0;
    VectorConstantCoefficient nx(nxVec), ny(nyVec), n1(n1Vec);

    MInv_->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));

    double alpha = -1.0, beta = 0.0;

    Kx_->AddDomainIntegrator(new DerivativeIntegrator(one, 0));
    Kx_->AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(nx, alpha, beta)));

    Ky_->AddDomainIntegrator(new DerivativeIntegrator(one, 1));
    Ky_->AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(ny, alpha, beta)));
}

void Solver::assembleBilinearForms()
{
    MInv_->Assemble();
    int skip_zeros = 0;
    Kx_->Assemble(skip_zeros);
    Ky_->Assemble(skip_zeros);
}

void Solver::finalizeBilinearForms()
{
    MInv_->Finalize();
    int skip_zeros = 0;
    Kx_->Finalize(skip_zeros);
    Ky_->Finalize(skip_zeros);
}


void Solver::run() {
    double t = 0.0;

    Vector aux(fes_->GetVSize());
    Vector ezNew(fes_->GetVSize());
    Vector hxNew(fes_->GetVSize());
    Vector hyNew(fes_->GetVSize());

    bool done = false;
    for (int ti = 0; !done; )
    {
        double dt_real = std::min(opts_.dt, opts_.t_final - t);

        // Update E.
        Kx_->Mult(hy_, aux);
        Ky_->AddMult(hx_, aux, -1.0);
        MInv_->Mult(aux, ezNew);
        ezNew *= -opts_.dt;
        ezNew.Add(1.0, ez_);


        // Update H.
        Kx_->Mult(ezNew, aux);
        MInv_->Mult(aux, hyNew);
        hyNew *= -opts_.dt;
        hyNew.Add(1.0, hy_);

        Ky_->Mult(ezNew, aux);
        MInv_->Mult(aux, hxNew);
        hxNew *= opts_.dt;
        hxNew.Add(1.0, hx_);

        ez_ = ezNew;
        hx_ = hxNew;
        hy_ = hyNew;

        t += opts_.dt;
        ti++;

        done = (t >= opts_.t_final - 1e-8 * opts_.dt);

        if (done || ti % opts_.vis_steps == 0)
        {
            std::cout << "time step: " << ti << ", time: " << t << std::endl;

            if (opts_.paraview)
            {
                pd_->SetCycle(ti);
                pd_->SetTime(t);
                pd_->Save();
            }
        }
    }
}
}
