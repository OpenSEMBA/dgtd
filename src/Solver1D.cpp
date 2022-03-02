#include "Solver1D.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell {

Solver1D::Solver1D(const Options& opts, const Mesh& mesh)
{
	checkOptionsAreValid(opts, mesh);

	mesh_ = mfem::Mesh(mesh, true);
	opts_ = opts;

	fec_ = std::make_unique<DG_FECollection>(opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	MInv_ = buildMassMatrix();
	Kx_ = buildDerivativeAndFluxOperator(X);

	Ez_.SetSpace(fes_.get());
	Ez_.ProjectCoefficient(ConstantCoefficient(0.0));

	Hy_.SetSpace(fes_.get());
	Hy_.ProjectCoefficient(ConstantCoefficient(0.0));

	initializeParaviewData();
}

void Solver1D::checkOptionsAreValid(const Options& opts, const Mesh& mesh)
{
	if (mesh.Dimension() != 1) {
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

std::unique_ptr<mfem::BilinearForm> Solver1D::buildMassMatrix() const
{
	auto MInv = std::make_unique<BilinearForm>(fes_.get());
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}

std::unique_ptr<mfem::BilinearForm> Solver1D::buildDerivativeAndFluxOperator(const Direction& d) const
{

	assert(d == X, "Incorrect argument for direction.");

	ConstantCoefficient one(1.0);

	auto K = std::make_unique<BilinearForm>(fes_.get());
	K->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, X)));

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	K->AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));
	K->AddBdrFaceIntegrator(
		new DGTraceIntegrator(n[0], -1.0, 0.0));

	K->Assemble();
	K->Finalize();

	return K;
}

void Solver1D::setInitialElectricField(std::function<ElectricField(const Position&)> f)
{
	Ez_.ProjectCoefficient(FunctionCoefficient(f));
}

void Solver1D::initializeParaviewData()
{
	pd_ = NULL;
	pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh_);
	pd_->SetPrefixPath("ParaView");
	pd_->RegisterField("Ez", &Ez_);
	pd_->RegisterField("Hy", &Hy_);
	pd_->SetLevelsOfDetail(opts_.order);
	pd_->SetDataFormat(VTKFormat::BINARY);
	opts_.order > 0 ? pd_->SetHighOrderOutput(true) : pd_->SetHighOrderOutput(false);
}

void Solver1D::run()
{
	double time = 0.0;
	bool done = false;

	Vector aux(fes_->GetVSize());
	Vector ezNew(fes_->GetVSize());
	Vector hyNew(fes_->GetVSize());

	pd_->SetCycle(0);
	pd_->SetTime(0.0);
	pd_->Save();

	for (int cycle = 0; !done;)
	{

		// Update E.
		Kx_->Mult(Hy_, aux);
		MInv_->Mult(aux, ezNew);
		ezNew *= opts_.dt;
		ezNew.Add(1.0, Ez_);

		// Update H.
		Kx_->Mult(Ez_, aux);
		MInv_->Mult(aux, hyNew);
		hyNew *= opts_.dt;
		hyNew.Add(1.0, Hy_);

		Ez_ = ezNew;
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