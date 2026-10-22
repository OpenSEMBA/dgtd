#pragma once

#include "Types.h"

#include <nlohmann/json.hpp>
#include <vector>

namespace maxwell {

/// Single-pole electric Lorentz on one mesh attribute.
/// Rates are SI values divided by c_SI (solver time is light-meters).
/// epsilon_ stored on Material for this tag is eps_inf.
struct LorentzProperties {
	Attribute geom_tag = 0;
	double eps_inf = 1.0;
	double omega_p = 0.0;
	double omega_1 = 0.0;
	double gamma = 0.0;
};

/// Parse a `lorentz` object. geom_tag is filled by the caller.
LorentzProperties parseLorentzObject(const nlohmann::json& lorentz_json);

} // namespace maxwell
