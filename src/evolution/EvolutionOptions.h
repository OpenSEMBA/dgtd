#pragma once

#include "components/Types.h"

namespace maxwell {

enum EvolutionOperatorType {
	Maxwell = 0,
	Global = 1,
	Hesthaven = 2
};

struct EvolutionOptions {
	EvolutionOperatorType op{ Global };
	int order{ 2 };
	double alpha = 1.0;
	bool spectral{ false };
};

}