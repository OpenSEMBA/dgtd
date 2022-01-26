#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

Solver::Solver() {
    initialize();
}

Solver::Solver(const Options& inOpts, const int nx, const int ny, mfem::Element::Type type, bool generateEdges) {
    opts = inOpts;
    initializeMesh();
    initializeFiniteElementSpace(initializeMesh());
    initializeBilinearForms(initializeFiniteElementSpace(initializeMesh()));
}

void Solver::initialize() {
    std::cout.precision(opts.precision);

    Device device(opts.device_config);

    Mesh mesh(opts.mesh_file, 1, 1);

    int dim = mesh.Dimension();

    for (int lev = 0; lev < opts.ref_levels; lev++)
    {
        mesh.UniformRefinement();
    }

    mesh.GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax, std::max(opts.order, 1));

    DG_FECollection fec(opts.order, dim, BasisType::GaussLobatto);
    fes = new FiniteElementSpace(&mesh, &fec);
    std::cout << "Number of unknowns per field:    " << fes->GetVSize() << std::endl;

    ConstantCoefficient zero(0.0), one(1.0), mOne(-1.0);
    Vector nxVec(2);  nxVec(0) = 1.0; nxVec(1) = 0.0;
    Vector nyVec(2);  nyVec(0) = 0.0; nyVec(1) = 1.0;
    Vector n1Vec(2);  n1Vec(0) = 1.0; n1Vec(1) = 1.0;
    VectorConstantCoefficient nx(nxVec), ny(nyVec), n1(n1Vec);

    MInv = new BilinearForm(fes);
    Kx = new BilinearForm(fes);
    Ky = new BilinearForm(fes);

    MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));

    double alpha = -1.0, beta = 0.0;

    Kx->AddDomainIntegrator(new DerivativeIntegrator(one, 0));
    Kx->AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(nx, alpha, beta)));

    Ky->AddDomainIntegrator(new DerivativeIntegrator(one, 1));
    Ky->AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(ny, alpha, beta)));

    MInv->Assemble();
    int skip_zeros = 0;
    Kx->Assemble(skip_zeros);
    Ky->Assemble(skip_zeros);

    MInv->Finalize();
    Kx->Finalize(skip_zeros);
    Ky->Finalize(skip_zeros);

    FunctionCoefficient u0(u0_function);
    GridFunction ez(fes), hx(fes), hy(fes);
    ez.ProjectCoefficient(u0);
    hx.ProjectCoefficient(zero);
    hy.ProjectCoefficient(zero);

    pd = NULL;
    if (opts.paraview)
    {
        pd = new ParaViewDataCollection("Example9", &mesh);
        pd->SetPrefixPath("ParaView");
        pd->RegisterField("ez", &ez);
        pd->RegisterField("hx", &hx);
        pd->RegisterField("hy", &hy);
        pd->SetLevelsOfDetail(opts.order);
        pd->SetDataFormat(VTKFormat::BINARY);
        opts.order > 0 ? pd->SetHighOrderOutput(true) : pd->SetHighOrderOutput(false);
        pd->SetCycle(0);
        pd->SetTime(0.0);
        pd->Save();
    }
}

FiniteElementSpace Solver::initializeFiniteElementSpace(Mesh mesh) {
    DG_FECollection fec(opts.order, mesh.Dimension(), BasisType::GaussLobatto);
    FiniteElementSpace fes(&mesh, &fec);
    std::cout << "Number of unknowns per field:    " << fes.GetTrueVSize() << std::endl;
    return fes;
}

void Solver::initializeBilinearForms(FiniteElementSpace fes) {

    ConstantCoefficient zero(0.0), one(1.0), mOne(-1.0);

    Vector nxVec(2);  nxVec(0) = 1.0; nxVec(1) = 0.0;
    Vector nyVec(2);  nyVec(0) = 0.0; nyVec(1) = 1.0;
    Vector n1Vec(2);  n1Vec(0) = 1.0; n1Vec(1) = 1.0;
    VectorConstantCoefficient nx(nxVec), ny(nyVec), n1(n1Vec);

    BilinearForm MInv(&fes);
    BilinearForm Kx(&fes);
    BilinearForm Ky(&fes);

    MInv.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));

    double alpha = -1.0, beta = 0.0;

    Kx.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
    Kx.AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(nx, alpha, beta)));

    Ky.AddDomainIntegrator(new DerivativeIntegrator(one, 1));
    Ky.AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(ny, alpha, beta)));

    MInv.Assemble();
    int skip_zeros = 0;
    Kx.Assemble(skip_zeros);
    Ky.Assemble(skip_zeros);

    MInv.Finalize();
    Kx.Finalize(skip_zeros);
    Ky.Finalize(skip_zeros);
}


void Solver::run() {
    double t = 0.0;

    Vector aux(fes->GetVSize());
    Vector ezNew(fes->GetVSize());
    Vector hxNew(fes->GetVSize());
    Vector hyNew(fes->GetVSize());

    bool done = false;
    for (int ti = 0; !done; )
    {
        double dt_real = std::min(opts.dt, opts.t_final - t);


        // Update E.
        Kx->Mult(hy, aux);
        Ky->AddMult(hx, aux, -1.0);
        MInv->Mult(aux, ezNew);
        ezNew *= -opts.dt;
        ezNew.Add(1.0, ez);


        // Update H.
        Kx->Mult(ezNew, aux);
        MInv->Mult(aux, hyNew);
        hyNew *= -opts.dt;
        hyNew.Add(1.0, hy);

        Ky->Mult(ezNew, aux);
        MInv->Mult(aux, hxNew);
        hxNew *= opts.dt;
        hxNew.Add(1.0, hx);

        ez = ezNew;
        hx = hxNew;
        hy = hyNew;

        t += opts.dt;
        ti++;

        done = (t >= opts.t_final - 1e-8 * opts.dt);

        if (done || ti % opts.vis_steps == 0)
        {
            std::cout << "time step: " << ti << ", time: " << t << std::endl;

            if (opts.paraview)
            {
                pd->SetCycle(ti);
                pd->SetTime(t);
                pd->Save();
            }
        }
    }

    delete pd;
}
