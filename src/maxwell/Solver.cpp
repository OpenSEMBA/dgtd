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
	probesManager_{ probes, fes_, fields_, opts_},
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

	checkOptionsAreValid(opts_);
}

void Solver::checkOptionsAreValid(const SolverOptions& opts)
{
	if ((opts.order < 0) ||
		(opts.t_final < 0)) {
		throw std::exception("Incorrect parameters in Options");
	}

	if (opts.dt == 0.0) {
		if (fes_.GetMesh()->Dimension() > 1) {
			throw std::exception("Automatic LTS calculation not implemented yet for Dimensions higher than 1.");
		}
		opts_.dt = calculateLTS();
	}

	for (const auto& bdrMarker : model_.getBoundaryToMarker())
	{
		if (bdrMarker.first == BdrCond::SMA && opts_.evolutionOperatorOptions.fluxType == FluxType::Centered) {
			throw std::exception("SMA and Centered FluxType are not compatible.");
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
	double res = std::numeric_limits<double>::max();
	for (int elemId = 0; elemId < fes.GetMesh()->ElementToElementTable().Size(); ++elemId) {
		Array<int> dofs;
		fes.GetElementDofs(elemId, dofs);
		for (int i = 0; i < dofs.Size(); ++i) {
			for (int j = i + 1; j < dofs.Size(); ++j) {
				res = std::min(res, std::abs(nodes[dofs[i]] - nodes[dofs[j]]));
			}
		}
	}
	return res;
}

double Solver::calculateLTS()
{
	double signalSpeed = 1.0;
	return (opts_.CFL * getMinimumInterNodeDistance(fes_)) / (pow(opts_.order, 1.5) * signalSpeed);
}

void Solver::run()
{
	while (time_ <= opts_.t_final - 1e-8*opts_.dt) {
		double truedt{ std::min(opts_.dt, opts_.t_final - time_) };
		odeSolver_->Step(fields_.allDOFs, time_, truedt);
		probesManager_.updateProbes(time_);
	}
}

}
