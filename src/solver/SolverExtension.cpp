#include "SolverExtension.h"

namespace maxwell {

struct SGBCContainer {



};

void SGBCSolver::initSGBCPreReqs()
{
	if (model_.getSGBCToMarker().find(BdrCond::SGBC) != model_.getSGBCToMarker().end()) {

		sourcesManager_.initTFSFPreReqs(model_.getConstMesh(), model_.getSGBCToMarker().at(BdrCond::SGBC));

		auto sub_fes = sourcesManager_.getGlobalTFSFSpace();
		auto sub_mesh = sub_fes->GetMesh();

		auto matrix{ assembleInteriorFluxMatrix(*sub_fes) };
		auto maps{ mapConnectivity(matrix.get()) };

		Mesh mesh1d = Mesh::MakeCartesian1D(10);
		DG_FECollection dgfec1d(3, 1, opts_.basisType);

		sgbcEvol_.resize(maps.first.size());
		

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