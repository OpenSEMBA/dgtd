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
	bool export_evolution_operator{ false };
	bool is_sgbc_solver = false;
	/// When true, GlobalEvolution builds only SGBC face ops (no vacuum A / PML / Debye /
	/// Lorentz / TFSF). Used as HesthavenEvolution's nested coating engine.
	bool hesthaven_sgbc_host = false;
};

}