#pragma once

#include "components/Types.h"

namespace maxwell {

struct EvolutionOptions {
	int order{ 2 };
	FluxType fluxType{ FluxType::Upwind };
	bool spectral{ false };
};

}