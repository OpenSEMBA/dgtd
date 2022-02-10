#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell {

Solver::Solver(const Options& opts, const Mesh& mesh) 
{

    checkOptionsAreValid(opts, mesh);

    mesh_ = mfem::Mesh(mesh, true);
    opts_ = opts;

    fec_ = std::make_unique<DG_FECollection>(opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
    fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
    fecH1_ = std::make_unique<H1_FECollection>(opts_.order, mesh_.Dimension());
    fesH1_ = std::make_unique<FiniteElementSpace>(&mesh_, fecH1_.get());

    boundaryTDOFs_ = buildEssentialTrueDOF();

    MInv_ = buildMassMatrix();
    Kx_ = buildDerivativeOperator(X);
    Ky_ = buildDerivativeOperator(Y);

    ez_.SetSpace(fes_.get());
    ez_.ProjectCoefficient(ConstantCoefficient(0.0));

    hx_.SetSpace(fes_.get());
    hx_.ProjectCoefficient(ConstantCoefficient(0.0));

    hy_.SetSpace(fes_.get());
    hy_.ProjectCoefficient(ConstantCoefficient(0.0));

    initializeParaviewData();

}

void Solver::checkOptionsAreValid(const Options& opts, const Mesh& mesh) 
{
    if (mesh.Dimension() != 2) {
        throw std::exception("Incorrect Dimension for mesh");
    }
    if ((opts.order < 0) ||
        (opts.t_final < 0) ||
        (opts.dt < 0) ||
        (opts.vis_steps < 1) ||
        (opts.precision < 1)) {
        throw std::exception("Incorrect parameters in Options");
    }
}

mfem::Array<int> Solver::buildEssentialTrueDOF()
{
    Array<int> ess_tdof_list;
    if (mesh_.bdr_attributes.Size())
    {
        Array<int> ess_bdr(mesh_.bdr_attributes.Max());
        ess_bdr = 1;
        fesH1_.get()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }
    std::cout << ess_tdof_list << std::endl;
    return ess_tdof_list;
}

std::unique_ptr<mfem::BilinearForm> Solver::buildMassMatrix() const
{
    auto MInv = std::make_unique<BilinearForm>(fes_.get());
    MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
    MInv->Assemble();
    MInv->Finalize();
    return MInv;
}

std::unique_ptr<mfem::BilinearForm> Solver::buildDerivativeOperator(const Direction& d) const
{
    assert(d != X || d != Y, "Incorrect argument for direction.");

    auto kDir = std::make_unique<BilinearForm>(fes_.get());

    ConstantCoefficient one(1.0);

    kDir->AddDomainIntegrator(
        new TransposeIntegrator(
            new DerivativeIntegrator(one, d)));

    std::vector<VectorConstantCoefficient> n = {
        VectorConstantCoefficient(Vector({1.0, 0.0})),
        VectorConstantCoefficient(Vector({0.0, 1.0}))
    };

    double alpha = 0.0;
    if (d == X) {
        alpha = 1.0;
    }
    else {
        alpha = -1.0;
    }

    double beta = 0.0;

    //kDir->AddBoundaryIntegrator(new DGTraceIntegrator(n[d], alpha, beta));

    //kDir->AddInteriorFaceIntegrator(
    //    new TransposeIntegrator(
    //        new DGTraceIntegrator(n[d], alpha, beta)));
    //kDir->AddBdrFaceIntegrator(
    //    new TransposeIntegrator(
    //        new DGTraceIntegrator(n[d], -alpha, beta)));


    int skip_zeros = 0;
    kDir->Assemble(skip_zeros);
    kDir->Finalize(skip_zeros);
    
    return kDir;
}

void Solver::setInitialElectricField(std::function<ElectricField(const Position&)> f) 
{
    ez_.ProjectCoefficient(FunctionCoefficient(f));
}

void Solver::initializeParaviewData()
{
    pd_ = NULL;
    pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView", &mesh_);
    pd_->SetPrefixPath("ParaView");
    pd_->RegisterField("ez", &ez_);
    pd_->RegisterField("hx", &hx_);
    pd_->RegisterField("hy", &hy_);
    pd_->SetLevelsOfDetail(opts_.order);
    pd_->SetDataFormat(VTKFormat::BINARY);
    opts_.order > 0 ? pd_->SetHighOrderOutput(true) : pd_->SetHighOrderOutput(false);
}


void Solver::run() 
{
    double time = 0.0;
    bool done = false;

    Vector aux(fes_->GetVSize());
    Vector ezNew(fes_->GetVSize());
    Vector hxNew(fes_->GetVSize());
    Vector hyNew(fes_->GetVSize());

    aux = 0.0; ezNew = 0.0; hxNew = 0.0; hyNew = 0.0;

    pd_->SetCycle(0);
    pd_->SetTime(0.0);
    pd_->Save();

    for (int cycle = 0; !done;)
    {

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

       //ez_.ProjectCoefficient(ConstantCoefficient(0.0), boundaryTDOFs_);

        time += opts_.dt;
        cycle++;

        done = (time >= opts_.t_final - 1e-8 * opts_.dt);

        if (done || cycle % opts_.vis_steps == 0) {
                pd_->SetCycle(cycle);
                pd_->SetTime(time);
                pd_->Save();
        }
    }
}

}
