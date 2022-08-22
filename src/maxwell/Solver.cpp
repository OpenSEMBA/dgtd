#include <fstream>
#include <iostream>
#include <algorithm>

#include "Solver.h"

using namespace mfem;

namespace maxwell {

FieldViews buildFieldsView(std::array<GridFunction, 3>& E, std::array<GridFunction, 3>& H)
{
	FieldViews r;
	for (const auto& x : { X, Y, Z }) {
		r.E[x] = &E[x];
		r.H[x] = &H[x];
	}
	return r;
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
	sources_{ sources },
	fec_{ opts_.order, model_.getMesh().Dimension(), BasisType::GaussLobatto},
	fes_{ &model_.getMesh(), &fec_ },
	maxwellEvol_{ &fes_, opts_.evolutionOperatorOptions, model_, sources_ }
{

	sol_ = Vector{
		FiniteElementEvolution::numberOfFieldComponents *
		FiniteElementEvolution::numberOfMaxDimensions *
		fes_.GetNDofs()
	};
	sol_ = 0.0;

	for (int d = X; d <= Z; d++) {
		E_[d].SetSpace(&fes_);
		E_[d].SetData(sol_.GetData() + d*fes_.GetNDofs());
		H_[d].SetSpace(&fes_);
		H_[d].SetData(sol_.GetData() + (d+3)*fes_.GetNDofs());
	}

	initializeFieldsFromSources();

	probesManager_ = ProbesManager{ probes, &fes_, buildFieldsView(E_, H_) };
}

void Solver::checkOptionsAreValid(const SolverOptions& opts)
{
	if ((opts.order < 0) ||
		(opts.t_final < 0) ||
		(opts.dt < 0)) {
		throw std::exception("Incorrect parameters in Options");
	}
}

void Solver::initializeFieldsFromSources()
{
	for (const auto& source: sources_) {
		std::function<double(const Position&)> f = 0;
		
		switch (model_.getConstMesh().Dimension()) {
		case 1:
			f = std::bind(&GaussianInitialField::evalGaussianFunction1D, &source, std::placeholders::_1);
			break;
		case 2:
			f = std::bind(&GaussianInitialField::evalGaussianFunction2D, &source, std::placeholders::_1);
			break;
		case 3:
			f = std::bind(&GaussianInitialField::evalGaussianFunction3D, &source, std::placeholders::_1);
			break;
		}

		switch (source.getFieldType()) {
		case FieldType::E:
			E_[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
			break;
		case FieldType::H:
			H_[source.getDirection()].ProjectCoefficient(FunctionCoefficient(f));
			break;
		}
	}
}

const GridFunction& Solver::getFieldInDirection(const FieldType& ft, const Direction& d) const
{
	assert(ft == E || ft == H);
	assert(d < 3);
	
	switch (ft) {
	case FieldType::E:
		return E_[d];
	case FieldType::H:
		return H_[d];
	}

	throw std::runtime_error("Invalid field type.");
}

void Solver::run()
{
	double time = 0.0;
	
	maxwellEvol_.SetTime(time);
	odeSolver_->Init(maxwellEvol_);
	
	probesManager_.updateProbes(time);
	
	bool done = false;
	while (!done) {
		odeSolver_->Step(sol_, time, opts_.dt);
		if (abs(time - opts_.t_final) < 1e-6) {
			done = true;
		}
		probesManager_.updateProbes(time);
	}

	probesManager_.updateProbes(time);
}

}
