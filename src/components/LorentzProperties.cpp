#include "LorentzProperties.h"

#include "math/PhysicalConstants.h"

#include <set>
#include <stdexcept>
#include <string>

namespace maxwell {

LorentzProperties parseLorentzObject(const nlohmann::json& lorentz_json)
{
	if (!lorentz_json.is_object()) {
		throw std::runtime_error(
			"lorentz must be an object with eps_inf, omega_p, omega_1, and gamma.");
	}

	const std::set<std::string> allowed{"eps_inf", "omega_p", "omega_1", "gamma"};
	for (auto it = lorentz_json.begin(); it != lorentz_json.end(); ++it) {
		if (allowed.count(it.key()) == 0) {
			throw std::runtime_error(
				"Unknown lorentz key '" + it.key() +
				"'. Expected eps_inf, omega_p, omega_1, and gamma.");
		}
	}
	for (const char* key : {"eps_inf", "omega_p", "omega_1", "gamma"}) {
		if (!lorentz_json.contains(key)) {
			throw std::runtime_error(std::string("lorentz is missing '") + key + "'.");
		}
	}

	LorentzProperties pole;
	pole.eps_inf = lorentz_json["eps_inf"].get<double>();
	const double omega_p_si = lorentz_json["omega_p"].get<double>();
	const double omega_1_si = lorentz_json["omega_1"].get<double>();
	const double gamma_si = lorentz_json["gamma"].get<double>();

	if (pole.eps_inf < 1.0) {
		throw std::runtime_error("lorentz eps_inf must be >= 1.");
	}
	if (!(omega_p_si > 0.0)) {
		throw std::runtime_error("lorentz omega_p must be > 0 (rad/s).");
	}
	if (omega_1_si < 0.0) {
		throw std::runtime_error("lorentz omega_1 must be >= 0 (rad/s).");
	}
	if (gamma_si < 0.0) {
		throw std::runtime_error("lorentz gamma must be >= 0 (rad/s).");
	}

	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;
	pole.omega_p = omega_p_si * inv_c;
	pole.omega_1 = omega_1_si * inv_c;
	pole.gamma = gamma_si * inv_c;
	return pole;
}

} // namespace maxwell
