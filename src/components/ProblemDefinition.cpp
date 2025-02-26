#include "ProblemDefinition.h"

namespace maxwell {
	ProblemDescription::ProblemDescription(Model& model, Probes& probes, Sources& sources, EvolutionOptions& opts) :
		model(model),
		probes(probes),
		sources(sources),
		opts(opts)
	{}

}