#pragma once

#include "components/Model.h"
#include "components/Probes.h"
#include "solver/SolverOptions.h"
#include "components/Sources.h"

namespace maxwell {

struct ProblemDescription {

	ProblemDescription(Model& model, Probes& probes, Sources& sources, EvolutionOptions& opts);

	Model model;
	Probes probes;
	Sources sources;
	EvolutionOptions opts;

};

}
