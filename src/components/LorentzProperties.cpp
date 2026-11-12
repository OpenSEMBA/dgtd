#include "LorentzProperties.h"

#include "math/PhysicalConstants.h"

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace maxwell {

namespace {

int parseAxisValue(const nlohmann::json& j, const std::string& context)
{
	if (j.is_number_integer()) {
		const int ax = j.get<int>();
		if (ax < 0 || ax > 2) {
			throw std::runtime_error(context + " axis must be 0, 1, or 2.");
		}
		return ax;
	}
	if (j.is_string()) {
		const std::string s = j.get<std::string>();
		if (s == "X" || s == "x") return 0;
		if (s == "Y" || s == "y") return 1;
		if (s == "Z" || s == "z") return 2;
	}
	throw std::runtime_error(context + " axis must be 0/1/2 or \"X\"/\"Y\"/\"Z\".");
}

Position parseCenterVector(const nlohmann::json& obj, const std::string& context)
{
	const std::string key = obj.contains("center") ? "center" : (obj.contains("mean") ? "mean" : "");
	if (key.empty()) {
		throw std::runtime_error(context + " profile requires 'center' (or 'mean').");
	}
	const auto& c = obj[key];
	if (c.is_number()) {
		Position pos(1);
		pos[0] = c.get<double>();
		return pos;
	}
	if (c.is_array() && !c.empty() && c.size() <= 3) {
		Position pos(static_cast<int>(c.size()));
		for (size_t i = 0; i < c.size(); ++i) {
			pos[static_cast<int>(i)] = c[i].get<double>();
		}
		return pos;
	}
	throw std::runtime_error(context + " 'center' must be a number or a 1..3 element array.");
}

struct ParsedRate {
	double scalar_solver = 0.0;
	std::shared_ptr<const Function> profile;
};

ParsedRate parseRateEntry(
	const nlohmann::json& entry,
	double unit_to_rad_s,
	bool require_positive,
	const std::string& label,
	const char* rate_unit)
{
	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;
	const double solver_scale = unit_to_rad_s * inv_c;

	if (entry.is_number()) {
		const double val_si = unit_to_rad_s * entry.get<double>();
		if (require_positive) {
			if (!(val_si > 0.0)) {
				throw std::runtime_error(label + " must be > 0 (" + rate_unit + ").");
			}
		} else {
			if (val_si < 0.0) {
				throw std::runtime_error(label + " must be >= 0 (" + rate_unit + ").");
			}
		}
		ParsedRate out;
		out.scalar_solver = val_si * inv_c;
		return out;
	}

	if (!entry.is_object()) {
		throw std::runtime_error(label + " must be a number or a profile object.");
	}
	if (!entry.contains("type") || !entry["type"].is_string()) {
		throw std::runtime_error(label + " profile object must define string 'type'.");
	}

	const std::string type = entry["type"].get<std::string>();
	ParsedRate out;

	if (type == "constant") {
		if (!entry.contains("value") || !entry["value"].is_number()) {
			throw std::runtime_error(label + " constant profile requires numeric 'value'.");
		}
		const double val_si = unit_to_rad_s * entry["value"].get<double>();
		if (require_positive && !(val_si > 0.0)) {
			throw std::runtime_error(label + " must be > 0 (" + rate_unit + ").");
		}
		if (!require_positive && val_si < 0.0) {
			throw std::runtime_error(label + " must be >= 0 (" + rate_unit + ").");
		}
		out.scalar_solver = val_si * inv_c;
		out.profile = std::make_shared<ConstantProfile>(out.scalar_solver);
		return out;
	}

	if (type == "gaussian" || type == "positive_gaussian_derivative") {
		const std::string peak_key = entry.contains("peak")
			? "peak"
			: (entry.contains("amplitude") ? "amplitude" : "");
		if (peak_key.empty() || !entry[peak_key].is_number()) {
			throw std::runtime_error(label + " " + type + " profile requires numeric 'peak' (or 'amplitude').");
		}
		const double peak_raw = entry[peak_key].get<double>();
		const double spread = entry.value("spread", 1.0);
		const double offset_raw = entry.value("offset", 0.0);
		const double cutoff_factor = entry.value("cutoff_factor", GaussianProfile::kDefaultCutoffFactor);
		if (require_positive && !(peak_raw > 0.0) && !(offset_raw > 0.0)) {
			throw std::runtime_error(label + " " + type + " profile requires positive peak or offset.");
		}
		if (peak_raw < 0.0 || offset_raw < 0.0) {
			throw std::runtime_error(label + " " + type + " profile 'peak' and 'offset' must be >= 0.");
		}
		if (!(spread > 0.0)) {
			throw std::runtime_error(label + " " + type + " profile 'spread' must be > 0.");
		}
		if (!(cutoff_factor > 0.0 && cutoff_factor < 1.0)) {
			throw std::runtime_error(label + " " + type + " profile 'cutoff_factor' must be in (0, 1).");
		}
		Position center(1);
		center[0] = 0.5;
		if (entry.contains("center") || entry.contains("mean")) {
			center = parseCenterVector(entry, label);
		}
		const double peak_solver = peak_raw * solver_scale;
		const double offset_solver = offset_raw * solver_scale;
		out.scalar_solver = offset_solver + peak_solver;
		if (type == "gaussian") {
			out.profile = std::make_shared<GaussianProfile>(
				peak_solver, spread, center, offset_solver, cutoff_factor);
		} else {
			const int axis = entry.contains("axis") ? parseAxisValue(entry["axis"], label) : 0;
			if (axis >= center.Size()) {
				throw std::runtime_error(label + " positive_gaussian_derivative axis exceeds center dimension.");
			}
			out.profile = std::make_shared<PositiveGaussianDerivativeProfile>(
				peak_solver, spread, center, axis, offset_solver, cutoff_factor);
		}
		return out;
	}

	if (type == "ramc_61km_omega_p" ||
	    type == "ramc_61km_nu_c_gupta" ||
	    type == "ramc_61km_nu_c_liu" ||
	    type == "ramc_50km_omega_p" ||
	    type == "ramc_50km_nu_c_gupta" ||
	    type == "ramc_50km_nu_c_liu") {
		const double nose_x = entry.contains("nose_x")
			? entry["nose_x"].get<double>()
			: entry.value("origin", 0.0);
		const int axis = entry.contains("axis") ? parseAxisValue(entry["axis"], label) : 0;
		const double direction = entry.value("direction", 1.0);
		if (direction != 1.0 && direction != -1.0) {
			throw std::runtime_error(label + " 'direction' must be 1.0 or -1.0.");
		}
		const bool ramc50 = type.rfind("ramc_50km_", 0) == 0;
		mfem::Vector probe_pt(1);
		probe_pt(0) = ramc50 ? -0.012 : -0.005;
		double ref_si = 0.0;
		if (!ramc50) {
			Ramc61kmPreset preset = Ramc61kmPreset::OmegaP;
			if (type == "ramc_61km_omega_p") {
				preset = Ramc61kmPreset::OmegaP;
				ref_si = ramc61kmStagnationOmegaP(probe_pt);
			} else if (type == "ramc_61km_nu_c_gupta") {
				preset = Ramc61kmPreset::NuC_Gupta;
				ref_si = ramc61kmStagnationNuC_Gupta(probe_pt);
			} else {
				preset = Ramc61kmPreset::NuC_Liu;
				ref_si = ramc61kmStagnationNuC_Liu(probe_pt);
			}
			out.profile = std::make_shared<Ramc61kmStagnationProfile>(
				preset, inv_c, nose_x, axis, direction);
		} else {
			Ramc50kmPreset preset = Ramc50kmPreset::OmegaP;
			if (type == "ramc_50km_omega_p") {
				preset = Ramc50kmPreset::OmegaP;
				ref_si = ramc50kmStagnationOmegaP(probe_pt);
			} else if (type == "ramc_50km_nu_c_gupta") {
				preset = Ramc50kmPreset::NuC_Gupta;
				ref_si = ramc50kmStagnationNuC_Gupta(probe_pt);
			} else {
				preset = Ramc50kmPreset::NuC_Liu;
				ref_si = ramc50kmStagnationNuC_Liu(probe_pt);
			}
			out.profile = std::make_shared<Ramc50kmStagnationProfile>(
				preset, inv_c, nose_x, axis, direction);
		}
		out.scalar_solver = ref_si * inv_c;
		return out;
	}

	throw std::runtime_error(
		"Unknown profile type '" + type + "' in " + label +
		". Supported: constant, gaussian, positive_gaussian_derivative, "
		"ramc_61km_omega_p, ramc_61km_nu_c_gupta, ramc_61km_nu_c_liu, "
		"ramc_50km_omega_p, ramc_50km_nu_c_gupta, ramc_50km_nu_c_liu.");
}

} // namespace

LorentzProperties parseLorentzObject(const nlohmann::json& lorentz_json)
{
	if (!lorentz_json.is_object()) {
		throw std::runtime_error(
			"lorentz must be an object with eps_inf and either "
			"(omega_p, omega_1, gamma/nu_c) in rad/s or (f_p, f_1, gamma/nu_c) in Hz.");
	}

	const std::set<std::string> allowed{
		"eps_inf", "omega_p", "omega_1", "gamma", "nu_c", "f_p", "f_1"};
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
	const bool has_gamma = lorentz_json.contains("gamma");
	const bool has_nu_c = lorentz_json.contains("nu_c");
	if (has_gamma && has_nu_c) {
		throw std::runtime_error("lorentz cannot define both 'gamma' and 'nu_c'.");
	}
	if (!has_gamma && !has_nu_c) {
		throw std::runtime_error("lorentz is missing 'gamma' (or 'nu_c').");
	}
	const std::string damp_key = has_gamma ? "gamma" : "nu_c";

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
	if (pole.eps_inf < 1.0) {
		throw std::runtime_error("lorentz eps_inf must be >= 1.");
	}

	if (hz_keys) {
		if (!has_f_p || !has_f_1) {
			throw std::runtime_error("lorentz Hz style requires f_p, f_1, and gamma.");
		}
		const double two_pi = 2.0 * M_PI;
		const ParsedRate wp = parseRateEntry(
			lorentz_json["f_p"], two_pi, true, "lorentz plasma frequency", "Hz");
		const ParsedRate w1 = parseRateEntry(
			lorentz_json["f_1"], two_pi, false, "lorentz resonance", "Hz");
		const ParsedRate gm = parseRateEntry(
			lorentz_json[damp_key], two_pi, false, "lorentz gamma", "Hz");
		pole.omega_p = wp.scalar_solver;
		pole.omega_1 = w1.scalar_solver;
		pole.gamma = gm.scalar_solver;
		pole.omega_p_profile = wp.profile;
		pole.omega_1_profile = w1.profile;
		pole.gamma_profile = gm.profile;
		return pole;
	}

	if (!has_omega_p || !has_omega_1) {
		throw std::runtime_error(
			"lorentz rad/s style requires omega_p, omega_1, and gamma.");
	}
	const ParsedRate wp = parseRateEntry(
		lorentz_json["omega_p"], 1.0, true, "lorentz plasma frequency", "rad/s");
	const ParsedRate w1 = parseRateEntry(
		lorentz_json["omega_1"], 1.0, false, "lorentz resonance", "rad/s");
	const ParsedRate gm = parseRateEntry(
		lorentz_json[damp_key], 1.0, false, "lorentz gamma", "rad/s");
	pole.omega_p = wp.scalar_solver;
	pole.omega_1 = w1.scalar_solver;
	pole.gamma = gm.scalar_solver;
	pole.omega_p_profile = wp.profile;
	pole.omega_1_profile = w1.profile;
	pole.gamma_profile = gm.profile;
	return pole;
}

LorentzProperties parseDrudeObject(const nlohmann::json& drude_json)
{
	if (!drude_json.is_object()) {
		throw std::runtime_error(
			"drude must be an object with optional eps_inf (default 1.0), "
			"(omega_p or f_p), and (nu_c or gamma).");
	}

	const std::set<std::string> allowed{
		"eps_inf", "omega_p", "f_p", "nu_c", "gamma"};
	for (auto it = drude_json.begin(); it != drude_json.end(); ++it) {
		if (allowed.count(it.key()) == 0) {
			throw std::runtime_error(
				"Unknown drude key '" + it.key() +
				"'. Use eps_inf with omega_p (rad/s) or f_p (Hz), and nu_c or gamma.");
		}
	}

	const bool has_f_p = drude_json.contains("f_p");
	const bool has_omega_p = drude_json.contains("omega_p");
	if (has_f_p && has_omega_p) {
		throw std::runtime_error("drude cannot mix 'f_p' (Hz) and 'omega_p' (rad/s).");
	}
	if (!has_f_p && !has_omega_p) {
		throw std::runtime_error("drude requires 'omega_p' (rad/s) or 'f_p' (Hz).");
	}

	const bool has_nu_c = drude_json.contains("nu_c");
	const bool has_gamma = drude_json.contains("gamma");
	if (has_nu_c && has_gamma) {
		throw std::runtime_error("drude cannot define both 'nu_c' and 'gamma'.");
	}
	if (!has_nu_c && !has_gamma) {
		throw std::runtime_error("drude requires 'nu_c' (or 'gamma').");
	}
	const std::string damp_key = has_nu_c ? "nu_c" : "gamma";

	LorentzProperties pole;
	pole.eps_inf = drude_json.value("eps_inf", 1.0);
	if (pole.eps_inf < 1.0) {
		throw std::runtime_error("drude eps_inf must be >= 1.");
	}

	const double unit_scale = has_f_p ? (2.0 * M_PI) : 1.0;
	const char* rate_unit = has_f_p ? "Hz" : "rad/s";
	const std::string wp_key = has_f_p ? "f_p" : "omega_p";

	const ParsedRate wp = parseRateEntry(
		drude_json[wp_key], unit_scale, true, "drude plasma frequency", rate_unit);
	const ParsedRate gm = parseRateEntry(
		drude_json[damp_key], unit_scale, false, "drude collision frequency", rate_unit);

	pole.omega_p = wp.scalar_solver;
	pole.omega_1 = 0.0;
	pole.gamma = gm.scalar_solver;
	pole.omega_p_profile = wp.profile;
	pole.omega_1_profile = nullptr;
	pole.gamma_profile = gm.profile;
	return pole;
}

} // namespace maxwell
