#include <fstream>
#include <iostream>
#include <algorithm>

#include "Solver.h"

using namespace mfem;

namespace maxwell {

std::unique_ptr<FiniteElementSpace> buildFiniteElementSpace(Mesh* m, FiniteElementCollection* fec)
{
	//if (dynamic_cast<ParMesh*>(m) != nullptr) {
	//	return std::make_unique<ParFiniteElementSpace>(m, fec);
	//}
	if (dynamic_cast<Mesh*>(m) != nullptr) {
		return std::make_unique<FiniteElementSpace>(m, fec);
	}
	throw std::runtime_error("Invalid mesh to build FiniteElementSpace");
}

Solver::Solver(const ProblemDescription& problem, const SolverOptions& options) :
	Solver(problem.model, problem.probes, problem.sources, options)
{}

Solver::Solver(
	const Model& model,
	const Probes& probes,
	const Sources& sources,
	const SolverOptions& options) :
	opts_{ options },
	model_{ model },
	fec_{ opts_.order, model_.getMesh().Dimension(), BasisType::GaussLobatto},
	fes_{ buildFiniteElementSpace(& model_.getMesh(), &fec_) },
	fields_{ *fes_ },
	sourcesManager_{ sources, *fes_ },
	probesManager_{ probes, *fes_, fields_, opts_ },
	time_{0.0}
{
	
	sourcesManager_.setInitialFields(fields_);
	switch (opts_.evolutionOperatorOptions.spectral) {
	case true:
		maxwellEvol_ = std::make_unique<MaxwellEvolution3D_Spectral>(
			*fes_, model_, sourcesManager_, opts_.evolutionOperatorOptions);
		break;
	default:
		maxwellEvol_ = std::make_unique<MaxwellEvolution3D>(
			*fes_, model_, sourcesManager_, opts_.evolutionOperatorOptions);
		break;
	}

	maxwellEvol_->SetTime(time_);
	odeSolver_->Init(*maxwellEvol_);

	probesManager_.updateProbes(time_);

	checkOptionsAreValid(opts_);

	if (opts_.dt == 0.0) {
		dt_ = getTimeStep();
	}
	else {
		dt_ = opts_.dt;
	}
}

void Solver::checkOptionsAreValid(const SolverOptions& opts) const
{
	if ((opts.order < 0) ||
		(opts.t_final < 0)) {
		throw std::runtime_error("Incorrect parameters in Options");
	}

	if (opts.dt == 0.0) {
		if (fes_->GetMesh()->Dimension() > 1) {
			throw std::runtime_error("Automatic TS calculation not implemented yet for Dimensions higher than 1.");
		}
	}

	for (const auto& bdrMarker : model_.getBoundaryToMarker())
	{
		if (bdrMarker.first == BdrCond::SMA && opts_.evolutionOperatorOptions.fluxType == FluxType::Centered) {
			throw std::runtime_error("SMA and Centered FluxType are not compatible.");
		}
	}
}

const PointProbe& Solver::getPointProbe(const std::size_t probe) const 
{ 
	return probesManager_.getPointProbe(probe); 
}

double getMinimumInterNodeDistance(FiniteElementSpace& fes)
{
	GridFunction nodes(&fes);
	fes.GetMesh()->GetNodes(nodes);
	double res{ std::numeric_limits<double>::max() };
	for (int e = 0; e < fes.GetMesh()->ElementToElementTable().Size(); ++e) {
		Array<int> dofs;
		fes.GetElementDofs(e, dofs);
		if (dofs.Size() == 1) {
			res = std::min(res, fes.GetMesh()->GetElementSize(e));
		}
		else {
			for (int i = 0; i < dofs.Size(); ++i) {
				for (int j = i + 1; j < dofs.Size(); ++j) {
					res = std::min(res, std::abs(nodes[dofs[i]] - nodes[dofs[j]]));
				}
			}
		}
	}
	return res;
}

double Solver::getTimeStep()
{
	double signalSpeed{ 1.0 };
	double maxTimeStep{ 0.0 };
	if (opts_.order == 0) {
		maxTimeStep = getMinimumInterNodeDistance(*fes_) / signalSpeed;
	}
	else {
		maxTimeStep = getMinimumInterNodeDistance(*fes_) / pow(opts_.order, 1.5) / signalSpeed;
	}
	return opts_.CFL * maxTimeStep;
}

void Solver::run()
{
	while (time_ <= opts_.t_final - 1e-8*dt_) {
		step();
	}
}

void Solver::step()
{
	double truedt{ std::min(dt_, opts_.t_final - time_) };
	odeSolver_->Step(fields_.allDOFs, time_, truedt);
	probesManager_.updateProbes(time_);
}

}
