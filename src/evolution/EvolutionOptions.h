#pragma once

#include "Types.h"

namespace maxwell {

struct EvolutionOptions {
	FluxType fluxType{ FluxType::Upwind };
	bool spectral{ false };
	bool eigenvals{ false };
	bool marketFile{ false };
	int powerMethod{ 0 };
};

}