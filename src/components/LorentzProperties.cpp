#include "LorentzProperties.h"

#include "math/PhysicalConstants.h"

#include <set>
#include <stdexcept>
#include <string>

namespace maxwell {

namespace {

void validateAndStore(
	LorentzProperties& pole,
	double omega_p_si,
	double omega_1_si,
	double gamma_si,
	const char* rate_unit)
{
	if (pole.eps_inf < 1.0) {
		throw std::runtime_error("lorentz eps_inf must be >= 1.");
	}
	if (!(omega_p_si > 0.0)) {
		throw std::runtime_error(
			std::string("lorentz plasma frequency must be > 0 (") + rate_unit + ").");
	}
	if (omega_1_si < 0.0) {
		throw std::runtime_error(
			std::string("lorentz resonance must be >= 0 (") + rate_unit + ").");
	}
	if (gamma_si < 0.0) {
		throw std::runtime_error(
			std::string("lorentz gamma must be >= 0 (") + rate_unit + ").");
	}

	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;
	pole.omega_p = omega_p_si * inv_c;
	pole.omega_1 = omega_1_si * inv_c;
	pole.gamma = gamma_si * inv_c;
}

} // namespace

LorentzProperties parseLorentzObject(const nlohmann::json& lorentz_json)
{
	if (!lorentz_json.is_object()) {
		throw std::runtime_error(
			"lorentz must be an object with eps_inf and either "
			"(omega_p, omega_1, gamma) in rad/s or (f_p, f_1, gamma) in Hz.");
	}

	const std::set<std::string> allowed{
		"eps_inf", "omega_p", "omega_1", "gamma", "f_p", "f_1"};
	for (auto it = lorentz_json.begin(); it != lorentz_json.end(); ++it) {
		if (allowed.count(it.key()) == 0) {
			throw std::runtime_error(
				"Unknown lorentz key '" + it.key() +
				"'. Use eps_inf with either omega_p/omega_1/gamma (rad/s) "
				"or f_p/f_1/gamma (Hz).");
		}
	}

	if (!lorentz_json.contains("eps_inf")) {
		throw std::runtime_error("lorentz is missing 'eps_inf'.");
	}
	if (!lorentz_json.contains("gamma")) {
		throw std::runtime_error("lorentz is missing 'gamma'.");
	}

	const bool has_f_p = lorentz_json.contains("f_p");
	const bool has_f_1 = lorentz_json.contains("f_1");
	const bool has_omega_p = lorentz_json.contains("omega_p");
	const bool has_omega_1 = lorentz_json.contains("omega_1");
	const bool hz_keys = has_f_p || has_f_1;
	const bool rad_keys = has_omega_p || has_omega_1;

	if (hz_keys && rad_keys) {
		throw std::runtime_error(
			"lorentz cannot mix Hz keys (f_p, f_1) with rad/s keys "
			"(omega_p, omega_1). Use one style; gamma takes that style's unit.");
	}

	LorentzProperties pole;
	pole.eps_inf = lorentz_json["eps_inf"].get<double>();

	if (hz_keys) {
		if (!has_f_p || !has_f_1) {
			throw std::runtime_error("lorentz Hz style requires f_p, f_1, and gamma.");
		}
		const double two_pi = 2.0 * M_PI;
		validateAndStore(
			pole,
			two_pi * lorentz_json["f_p"].get<double>(),
			two_pi * lorentz_json["f_1"].get<double>(),
			two_pi * lorentz_json["gamma"].get<double>(),
			"Hz");
		return pole;
	}

	if (!has_omega_p || !has_omega_1) {
		throw std::runtime_error(
			"lorentz rad/s style requires omega_p, omega_1, and gamma.");
	}
	validateAndStore(
		pole,
		lorentz_json["omega_p"].get<double>(),
		lorentz_json["omega_1"].get<double>(),
		lorentz_json["gamma"].get<double>(),
		"rad/s");
	return pole;
}

} // namespace maxwell
