#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell1D {

Solver::Solver(const Options& opts, const Mesh& mesh)
{
	checkOptionsAreValid(opts, mesh);

	mesh_ = mfem::Mesh(mesh, true);
	opts_ = opts;

	fec_ = std::make_unique<DG_FECollection>(opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	MInv_ = buildInverseMassMatrix();
	KxE_ = buildDerivativeAndFluxOperator(X, Electric);
	KxH_ = buildDerivativeAndFluxOperator(X, Magnetic);

	odeSolverE_ = new RK4Solver;
	odeSolverH_ = new RK4Solver;
		
	Ez_.SetSpace(fes_.get());
	Ez_.ProjectCoefficient(ConstantCoefficient(0.0));

	Hy_.SetSpace(fes_.get());
	Hy_.ProjectCoefficient(ConstantCoefficient(0.0));

	initializeParaviewData();
}

void Solver::checkOptionsAreValid(const Options& opts, const Mesh& mesh)
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

mfem::Array<int> Solver::buildEssentialTrueDOF()
{
	Array<int> ess_tdof_list;
	if (mesh_.bdr_attributes.Size())
	{
		Array<int> ess_bdr(mesh_.bdr_attributes.Max());
		ess_bdr = 1;
		fes_.get()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
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

std::unique_ptr<mfem::BilinearForm> Solver::buildDerivativeAndFluxOperator(
	const Direction& d, const FieldType& ft) const
{
	assert(d == X, "Incorrect argument for direction.");

	auto K = std::make_unique<BilinearForm>(fes_.get());
		
	ConstantCoefficient one(1.0);

	K->AddDomainIntegrator(
		new TransposeIntegrator(
			new DerivativeIntegrator(one, d)));

	std::vector<VectorConstantCoefficient> n = {
		VectorConstantCoefficient(Vector({1.0})),
	};

	double alpha;
	double beta;

	if (ft == Electric)
	{
		alpha = -1.0;
		beta = 0.0;
		K->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		K->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}
	else
	{
		alpha = -1.0;
		beta = 0.0;
		K->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
		K->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[d], alpha, beta));
	}

	K->Assemble();
	K->Finalize();

	return K;
}

void Solver::setInitialElectricField(std::function<ElectricField(const Position&)> f)
{
	Ez_.ProjectCoefficient(FunctionCoefficient(f));
}

void Solver::initializeParaviewData()
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

void Solver::run()
{
	FE_Evolution electricOp(MInv_, KxH_);
	FE_Evolution magneticOp(MInv_, KxE_);

	double t = 0.0;
	electricOp.SetTime(t);
	magneticOp.SetTime(t);
	odeSolverE_->Init(electricOp);
	odeSolverH_->Init(magneticOp);

	pd_->SetCycle(0);
	pd_->SetTime(0.0);
	pd_->Save();

	double time = 0.0;
	bool done = false;
	for (int cycle = 0; !done;)
	{
		double dt_real = std::min(opts_.dt, opts_.t_final - t);

		odeSolverE_->Step(Hy_, t, dt_real);
		odeSolverH_->Step(Ez_, t, dt_real);

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