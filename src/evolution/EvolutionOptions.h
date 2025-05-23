#pragma once

#include "components/Types.h"

namespace maxwell {

struct EvolutionOptions {
	int order{ 2 };
	FluxType fluxType{ FluxType::Upwind };
	double alpha = 1.0;
	bool spectral{ false };
};

}