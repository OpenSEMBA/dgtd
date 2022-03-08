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

    MInv_ = buildInverseMassMatrix();
    KxE_ = buildDerivativeAndFluxOperator(X,Electric);
    KyE_ = buildDerivativeAndFluxOperator(Y,Electric);
    KxH_ = buildDerivativeAndFluxOperator(X,Magnetic);
    KyH_ = buildDerivativeAndFluxOperator(Y,Magnetic);

    Ez_.SetSpace(fes_.get());
    Ez_.ProjectCoefficient(ConstantCoefficient(0.0));

    Hx_.SetSpace(fes_.get());
    Hx_.ProjectCoefficient(ConstantCoefficient(0.0));

    Hy_.SetSpace(fes_.get());
    Hy_.ProjectCoefficient(ConstantCoefficient(0.0));

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
        fes_.get()->GetBoundaryTrueDofs(ess_tdof_list);
    }
    return ess_tdof_list;
}

std::unique_ptr<mfem::BilinearForm> Solver::buildInverseMassMatrix() const
{
    auto MInv = std::make_unique<BilinearForm>(fes_.get());
    MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
    MInv->Assemble();
    MInv->Finalize();
    return MInv;
}

std::unique_ptr<mfem::BilinearForm> Solver::buildDerivativeAndFluxOperator(const Direction& d, const FieldType& ft) const
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

    double alpha;
    double beta;
    
    if (ft == Electric) 
    {
        alpha = -1.0;
        beta = 0.0;
        kDir->AddInteriorFaceIntegrator(
            new DGTraceIntegrator(n[d], alpha, beta));
        kDir->AddBdrFaceIntegrator(
            new DGTraceIntegrator(n[d], 0.0*alpha, beta));
    }
    else
    {
        alpha = -1.0;
        beta = 0.0;
        kDir->AddInteriorFaceIntegrator(
            new DGTraceIntegrator(n[d], alpha, beta));
        kDir->AddBdrFaceIntegrator(
            new DGTraceIntegrator(n[d], 2.0*alpha, beta));
    }

    int skip_zeros = 0;
    kDir->Assemble(skip_zeros);
    kDir->Finalize(skip_zeros);
    
    return kDir;
}

void Solver::setInitialElectricField(std::function<ElectricField(const Position&)> f) 
{
    Ez_.ProjectCoefficient(FunctionCoefficient(f));
}

void Solver::initializeParaviewData()
{
    pd_ = NULL;
    pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView", &mesh_);
    pd_->SetPrefixPath("ParaView");
    pd_->RegisterField("ez", &Ez_);
    pd_->RegisterField("hx", &Hx_);
    pd_->RegisterField("hy", &Hy_);
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

    pd_->SetCycle(0);
    pd_->SetTime(0.0);
    pd_->Save();

    for (int cycle = 0; !done;)
    {

        // Update E.
        KxH_->Mult(Hy_, aux);
        //Ky_->Mult(Hx_, aux);
        MInv_->Mult(aux, ezNew);
        ezNew *= -opts_.dt;
        ezNew.Add(1.0, Ez_);

        // Update H.
        KxE_->Mult(ezNew, aux);
        MInv_->Mult(aux, hyNew);
        hyNew *= -opts_.dt;
        hyNew.Add(1.0, Hy_);

        //Ky_->Mult(ezNew, aux);
        //MInv_->Mult(aux, hxNew);
        //hxNew *= opts_.dt;
        //hxNew.Add(1.0, Hx_);

        Ez_ = ezNew;
        /*Hx_ = hxNew;*/
        Hy_ = hyNew;

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
