#include "Solver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace maxwell1D {

Solver::Solver(const Options& opts, const Mesh& mesh)
{
	checkOptionsAreValid(opts, mesh);

	mesh_ = mfem::Mesh(mesh, true);
	opts_ = opts;

	fec_ = std::make_unique<DG_FECollection>(
		opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
	fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

	odeSolver_ = std::make_unique<RK4Solver>();

	maxwellEvol_ = std::make_unique<FE_Evolution>(fes_.get(), opts.evolutionOperatorOptions);
	
	sol_ = Vector(FE_Evolution::numberOfFieldComponents * fes_->GetNDofs());
	sol_ = 0.0;

	E_.SetSpace(fes_.get());
	E_.SetData(sol_.GetData());
	H_.SetSpace(fes_.get());
	H_.SetData(sol_.GetData() + fes_->GetNDofs());

	if (opts_.paraview) {
		initializeParaviewData();
	}
	if (opts_.glvis) {
		initializeGLVISData();
	}
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

void Solver::setInitialField(const FieldType& ft, std::function<double(const Position&)> f)
{
	switch (ft) {
	case FieldType::Electric:
		E_.ProjectCoefficient(FunctionCoefficient(f));
		return;
	case FieldType::Magnetic:
		H_.ProjectCoefficient(FunctionCoefficient(f));
		return;
	}	
}

const GridFunction& Solver::getField(const FieldType& ft) const
{
	switch (ft) {
	case FieldType::Electric:
		return E_;
	case FieldType::Magnetic:
		return H_;
	}
}

void Solver::initializeParaviewData()
{
	pd_ = NULL;
	pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh_);
	pd_->SetPrefixPath("ParaView");
	pd_->RegisterField("E", &E_);
	pd_->RegisterField("H", &H_);
	pd_->SetLevelsOfDetail(opts_.order);
	pd_->SetDataFormat(VTKFormat::BINARY);
	opts_.order > 0 ? pd_->SetHighOrderOutput(true) : pd_->SetHighOrderOutput(false);
}

void Solver::initializeGLVISData() 
{
	char vishost[] = "localhost";
	int  visport = 19916;
	sout_.open(vishost, visport);
	sout_.precision(opts_.precision);
	sout_ << "solution\n" << mesh_ << E_;
	sout_ << "pause\n";
	sout_ << std::flush;
	std::cout << "GLVis visualization paused."
		<< " Press space (in the GLVis window) to resume it.\n";
}

void Solver::storeInitialVisualizationValues() 
{
	if (opts_.paraview) {
		pd_->SetCycle(0);
		pd_->SetTime(0.0);
		pd_->Save();
	}
	if (opts_.glvis) {
		std::ofstream omesh("Maxwell1D_RK4.mesh");
		omesh.precision(opts_.precision);
		mesh_.Print(omesh);
		std::ofstream eSol("Maxwell1D_RK4-init-E.gf");
		eSol.precision(opts_.precision);
		E_.Save(eSol);
		std::ofstream hSol("Maxwell1D_RK4-init-H.gf");
		hSol.precision(opts_.precision);
		H_.Save(hSol);
	}
}

void Solver::run()
{
	Vector vals(2 * fes_.get()->GetVSize());
	
	double time = 0.0;

	maxwellEvol_->SetTime(time);
	odeSolver_->Init(*maxwellEvol_);

	storeInitialVisualizationValues();

	bool done = false;
	int cycle = 0;

	while (!done) {
		odeSolver_->Step(sol_, time, opts_.dt);

		done = (time >= opts_.t_final);

		cycle++;

		if (done || cycle % opts_.vis_steps == 0) {
			if (opts_.paraview) {
				pd_->SetCycle(cycle);
				pd_->SetTime(time);
				pd_->Save();
			}
			if (opts_.glvis) {
				sout_ << "solution\n" << mesh_ << E_ << std::flush;
			}
		}
	}
}
}