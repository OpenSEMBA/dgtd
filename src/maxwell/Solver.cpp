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
	time_{0.0},
	maxwellEvol_{ &fes_, opts_.evolutionOperatorOptions, model_, sourcesManager_.sources }
{
	sourcesManager_.setFields(fields_);

	maxwellEvol_.SetTime(time_);
	odeSolver_->Init(maxwellEvol_);

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

const PointsProbe& Solver::getPointsProbe(const std::size_t probe) const 
{ 
	return probesManager_.getPointsProbe(probe); 
}

const GridFunction& Solver::getFieldInDirection(const FieldType& ft, const Direction& d) const
{
	assert(ft == E || ft == H);
	assert(d < 3);
	
	switch (ft) {
	case FieldType::E:
		return fields_.E[d];
	case FieldType::H:
		return fields_.H[d];
	}

	throw std::runtime_error("Invalid field type.");
}


void Solver::run()
{
	while (time_ < opts_.t_final) {
		odeSolver_->Step(fields_.allDOFs, time_, opts_.dt);
		probesManager_.updateProbes(time_);
	}
}

}
