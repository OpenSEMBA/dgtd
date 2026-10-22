#pragma once

#include "Types.h"

#include <nlohmann/json.hpp>
#include <vector>

namespace maxwell {

/// Single-pole electric Debye on one mesh attribute.
/// tau_solver is τ_SI * c_SI (solver time is light-meters). epsilon_ stored on
/// Material for this tag is eps_inf, the instantaneous mass weight.
struct DebyeProperties {
	Attribute geom_tag = 0;
	double eps_inf = 1.0;
	double eps_s = 1.0;
	double tau_solver = 0.0;
};

/// Parse a `debye` object. geom_tag is filled by the caller.
DebyeProperties parseDebyeObject(const nlohmann::json& debye_json);

} // namespace maxwell
