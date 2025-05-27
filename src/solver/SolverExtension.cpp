#include "SolverExtension.h"

namespace maxwell {

struct SGBCContainer {



};

std::vector<NodePair> findDofPairs(BilinearForm& form) {
	std::vector<NodePair> npvec;
	double tol = 1e-8;
	for (auto r = 0; r < form.NumRows(); r++) {
		NodePair np;
		for (auto c = 0; c < form.NumCols(); c++) {
			const auto& val = form.Elem(r, c);
			if (std::abs(form.Elem(r, c)) >= tol && c != form.NumCols()) {
				for (auto c2 = c; c2 < form.NumCols(); c2++) {
					if (val == form.Elem(r, c2)) {
						np.first = r;
						np.second = c2;
						npvec.push_back(np);
						c2 = form.NumCols();
						c = form.NumCols();
					}
				}
			}
		}
	}
	return npvec;
}

void SGBCSolver::initSGBCPreReqs()
{
	if (model_.getSGBCToMarker().find(BdrCond::SGBC) != model_.getSGBCToMarker().end()) {

		BilinearForm sgbcform(fes_.get());
		sgbcform.AddInternalBoundaryFaceIntegrator(new mfemExtension::MaxwellDGInteriorJumpIntegrator({}, 1.0), model_.getSGBCToMarker().at(BdrCond::SGBC));
		sgbcform.Assemble();
		sgbcform.Finalize();

		std::vector<NodePair> npvec = findDofPairs(sgbcform);
		

	}
}

SGBCSolver::SGBCSolver(
	const Model& model,
	const Probes& probes,
	const Sources& sources,
	const SolverOptions& options) :
	opts_{ options },
	model_{ model },
	fec_{ opts_.evolution.order, model_.getMesh().Dimension(), opts_.basisType },
	fes_{ buildFiniteElementSpace(&model_.getMesh(), &fec_) },
	fields_{ *fes_ },
	sourcesManager_{ sources, *fes_, fields_ },
	probesManager_{ probes , *fes_, fields_, opts_ },
	time_{ 0.0 }
{

	checkOptionsAreValid(opts_);

	if (opts_.evolution.spectral == true) {
		performSpectralAnalysis(*fes_.get(), model_, opts_.evolution);
	}

	initSGBCPreReqs();

	maxwellEvol_ = assignEvolutionOperator();
	maxwellEvol_->SetTime(time_);

	if (opts_.timeStep == 0.0) {
		dt_ = estimateTimeStep(model_, opts_, *fes_, maxwellEvol_.get());
	}
	else {
		dt_ = opts_.timeStep;
	}

	odeSolver_->Init(*maxwellEvol_);

	probesManager_.updateProbes(time_);

}

}