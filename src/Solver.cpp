#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell {

Solver::Solver(const Options& opts, Mesh& mesh) 
{
    opts_ = opts;

    Device device(opts.device_config);
    //mesh.getBoundingBox...

    initializeFiniteElementSpace(mesh);

    initializeBilinearForms();
}

void Solver::initializeFiniteElementSpace(Mesh& mesh) 
{
    DG_FECollection fec(opts_.order, mesh.Dimension(), BasisType::GaussLobatto);
    fes_ = std::make_unique<FiniteElementSpace>(&mesh, &fec);
}

void Solver::initializeBilinearForms() 
{

    ConstantCoefficient zero(0.0), one(1.0), mOne(-1.0);

    Vector nxVec(2);  nxVec(0) = 1.0; nxVec(1) = 0.0;
    Vector nyVec(2);  nyVec(0) = 0.0; nyVec(1) = 1.0;
    Vector n1Vec(2);  n1Vec(0) = 1.0; n1Vec(1) = 1.0;
    VectorConstantCoefficient nx(nxVec), ny(nyVec), n1(n1Vec);

    MInv_ = std::make_unique<BilinearForm>(fes_.get());
    //BilinearForm Kx(fes_);
    //BilinearForm Ky(fes_);

    //MInv.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));

    //double alpha = -1.0, beta = 0.0;

    //Kx.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
    //Kx.AddInteriorFaceIntegrator(
    //    new TransposeIntegrator(new DGTraceIntegrator(nx, alpha, beta)));

    //Ky.AddDomainIntegrator(new DerivativeIntegrator(one, 1));
    //Ky.AddInteriorFaceIntegrator(
    //    new TransposeIntegrator(new DGTraceIntegrator(ny, alpha, beta)));

    //MInv.Assemble();
    //int skip_zeros = 0;
    //Kx.Assemble(skip_zeros);
    //Ky.Assemble(skip_zeros);

    //MInv.Finalize();
    //Kx.Finalize(skip_zeros);
    //Ky.Finalize(skip_zeros);
}


//void Solver::run() {
//    double t = 0.0;
//
//    Vector aux(fes->GetVSize());
//    Vector ezNew(fes->GetVSize());
//    Vector hxNew(fes->GetVSize());
//    Vector hyNew(fes->GetVSize());
//
//    bool done = false;
//    for (int ti = 0; !done; )
//    {
//        double dt_real = std::min(opts.dt, opts.t_final - t);
//
//
//        // Update E.
//        Kx->Mult(hy, aux);
//        Ky->AddMult(hx, aux, -1.0);
//        MInv->Mult(aux, ezNew);
//        ezNew *= -opts.dt;
//        ezNew.Add(1.0, ez);
//
//
//        // Update H.
//        Kx->Mult(ezNew, aux);
//        MInv->Mult(aux, hyNew);
//        hyNew *= -opts.dt;
//        hyNew.Add(1.0, hy);
//
//        Ky->Mult(ezNew, aux);
//        MInv->Mult(aux, hxNew);
//        hxNew *= opts.dt;
//        hxNew.Add(1.0, hx);
//
//        ez = ezNew;
//        hx = hxNew;
//        hy = hyNew;
//
//        t += opts.dt;
//        ti++;
//
//        done = (t >= opts.t_final - 1e-8 * opts.dt);
//
//        if (done || ti % opts.vis_steps == 0)
//        {
//            std::cout << "time step: " << ti << ", time: " << t << std::endl;
//
//            if (opts.paraview)
//            {
//                pd->SetCycle(ti);
//                pd->SetTime(t);
//                pd->Save();
//            }
//        }
//    }
//
//    delete pd;
//}
}