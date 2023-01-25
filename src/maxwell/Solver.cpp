#include <fstream>
#include <iostream>
#include <algorithm>

#include "Solver.h"

using namespace mfem;

namespace maxwell {

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
	fes_{ &model_.getMesh(), &fec_ },
	fields_{ fes_ },
	sourcesManager_{ sources, fes_ },
	probesManager_{ probes, fes_, fields_},
	time_{0.0}
{
	sourcesManager_.setInitialFields(fields_);
	switch (fes_.GetMesh()->Dimension()) {
	case 1:
		maxwellEvol_ = std::make_unique<MaxwellEvolution1D>(fes_, model_, sourcesManager_, opts_.evolutionOperatorOptions);
		break;
	case 2:
		maxwellEvol_ = std::make_unique<MaxwellEvolution2D>(fes_, model_, sourcesManager_, opts_.evolutionOperatorOptions);
		break;
	default:
		maxwellEvol_ = std::make_unique<MaxwellEvolution3D>(fes_, model_, sourcesManager_, opts_.evolutionOperatorOptions);
		break;
	}
	maxwellEvol_->SetTime(time_);
	odeSolver_->Init(*maxwellEvol_);

	probesManager_.updateProbes(time_);
}

void Solver::checkOptionsAreValid(const SolverOptions& opts)
{
	if ((opts.order < 0) ||
		(opts.t_final < 0) ||
		(opts.dt < 0)) {
		throw std::exception("Incorrect parameters in Options");
	}
}

const PointProbe& Solver::getPointProbe(const std::size_t probe) const 
{ 
	return probesManager_.getPointProbe(probe); 
}

void Solver::run()
{
	while ( std::abs(time_ - opts_.t_final) < 1e-6 || time_ < opts_.t_final) {
		odeSolver_->Step(fields_.allDOFs, time_, opts_.dt);
		probesManager_.updateProbes(time_);
	}
}

}
