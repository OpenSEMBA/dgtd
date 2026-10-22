#include "DebyeProperties.h"

#include "math/PhysicalConstants.h"

#include <set>
#include <stdexcept>
#include <string>

namespace maxwell {

DebyeProperties parseDebyeObject(const nlohmann::json& debye_json)
{
	if (!debye_json.is_object()) {
		throw std::runtime_error("debye must be an object with eps_inf, eps_s, and tau.");
	}

	const std::set<std::string> allowed{"eps_inf", "eps_s", "tau"};
	for (auto it = debye_json.begin(); it != debye_json.end(); ++it) {
		if (allowed.count(it.key()) == 0) {
			throw std::runtime_error(
				"Unknown debye key '" + it.key() + "'. Expected eps_inf, eps_s, and tau.");
		}
	}
	for (const char* key : {"eps_inf", "eps_s", "tau"}) {
		if (!debye_json.contains(key)) {
			throw std::runtime_error(
				std::string("debye is missing '") + key + "'.");
		}
	}

	DebyeProperties pole;
	pole.eps_inf = debye_json["eps_inf"].get<double>();
	pole.eps_s = debye_json["eps_s"].get<double>();
	const double tau_si = debye_json["tau"].get<double>();

	if (pole.eps_inf < 1.0) {
		throw std::runtime_error("debye eps_inf must be >= 1.");
	}
	if (!(pole.eps_s > pole.eps_inf)) {
		throw std::runtime_error("debye eps_s must be greater than eps_inf.");
	}
	if (!(tau_si > 0.0)) {
		throw std::runtime_error("debye tau must be > 0 (seconds).");
	}

	pole.tau_solver = tau_si * physicalConstants::speedOfLight_SI;
	return pole;
}

} // namespace maxwell
