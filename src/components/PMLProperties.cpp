#include "PMLProperties.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>

namespace maxwell {

namespace {

Direction parseAxisToken(const std::string& token)
{
	if (token == "X" || token == "x") {
		return X;
	}
	if (token == "Y" || token == "y") {
		return Y;
	}
	if (token == "Z" || token == "z") {
		return Z;
	}
	throw std::runtime_error("PML active_axes entry must be \"X\", \"Y\", or \"Z\". Got: " + token);
}

PMLStretchMode parseStretchMode(const nlohmann::json& mat_json)
{
	if (!mat_json.contains("stretch_mode")) {
		return PMLStretchMode::Box;
	}
	const auto& v = mat_json["stretch_mode"];
	if (v.is_number_integer()) {
		const int m = v.get<int>();
		if (m == 0) {
			return PMLStretchMode::Box;
		}
		if (m == 1) {
			return PMLStretchMode::Radial;
		}
		throw std::runtime_error("PML stretch_mode integer must be 0 (box) or 1 (radial).");
	}
	if (v.is_string()) {
		std::string s = v.get<std::string>();
		std::transform(s.begin(), s.end(), s.begin(),
		               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		if (s == "box") {
			return PMLStretchMode::Box;
		}
		if (s == "radial") {
			return PMLStretchMode::Radial;
		}
		throw std::runtime_error(
			"PML stretch_mode must be \"box\" or \"radial\" (or 0/1). Got: " +
			v.get<std::string>());
	}
	throw std::runtime_error("PML stretch_mode must be a string or integer.");
}

std::optional<std::array<double, 3>> parseRadialCenter(const nlohmann::json& mat_json,
                                                       int mesh_dim)
{
	if (!mat_json.contains("radial_center")) {
		return std::nullopt;
	}
	const auto& arr = mat_json["radial_center"];
	if (!arr.is_array() || static_cast<int>(arr.size()) < mesh_dim ||
	    static_cast<int>(arr.size()) > 3) {
		throw std::runtime_error(
			"PML radial_center must be an array of length mesh_dim..3.");
	}
	std::array<double, 3> c{{0.0, 0.0, 0.0}};
	for (size_t i = 0; i < arr.size(); ++i) {
		c[i] = arr[i].get<double>();
	}
	return c;
}

} // namespace

std::set<Direction> parseActiveAxes(const nlohmann::json& mat_json, int mesh_dim)
{
	std::set<Direction> axes;
	if (!mat_json.contains("active_axes")) {
		throw std::runtime_error("PML material block requires 'active_axes'.");
	}
	for (const auto& entry : mat_json["active_axes"]) {
		const Direction d = parseAxisToken(entry.get<std::string>());
		if (d >= mesh_dim) {
			throw std::runtime_error(
				"PML active_axes direction exceeds mesh dimension.");
		}
		axes.insert(d);
	}
	if (axes.empty()) {
		throw std::runtime_error("PML active_axes must contain at least one direction.");
	}
	return axes;
}

void validatePMLMaterialBlock(const nlohmann::json& mat_json)
{
	if (mat_json.contains("bulk_conductivity")) {
		throw std::runtime_error(
			"PML material must not define bulk_conductivity. Use stretch profiles instead.");
	}
	if (mat_json.contains("relative_permittivity") ||
	    mat_json.contains("relative_permeability")) {
		throw std::runtime_error(
			"PML material must not define relative_permittivity/permeability. "
			"Use matches_vacuum instead.");
	}
	if (mat_json.contains("matches_vacuum") && !mat_json["matches_vacuum"].get<bool>()) {
		throw std::runtime_error("Only matches_vacuum: true is supported for volumetric PML.");
	}
	if (mat_json.contains("alpha_max") && mat_json["alpha_max"].get<double>() != 0.0) {
		throw std::runtime_error(
			"PML alpha_max > 0 is not supported yet (CFS pole deferred). Use alpha_max: 0 "
			"or omit the field.");
	}
}

PMLProperties parsePMLMaterialBlock(const nlohmann::json& mat_json, int mesh_dim)
{
	validatePMLMaterialBlock(mat_json);

	PMLProperties props;
	props.matches_vacuum = mat_json.value("matches_vacuum", true);
	props.grading_order = mat_json.value("grading_order", 3);
	props.target_reflection = mat_json.value("target_reflection", 1e-6);
	props.active_axes = parseActiveAxes(mat_json, mesh_dim);
	props.stretch_mode = parseStretchMode(mat_json);
	props.radial_center = parseRadialCenter(mat_json, mesh_dim);
	props.kappa_max = mat_json.value("kappa_max", 1.0);
	props.alpha_max = mat_json.value("alpha_max", 0.0);

	if (props.kappa_max < 1.0) {
		throw std::runtime_error("PML kappa_max must be >= 1.");
	}
	if (props.alpha_max < 0.0) {
		throw std::runtime_error("PML alpha_max must be >= 0.");
	}

	if (props.grading_order < 0) {
		throw std::runtime_error("PML grading_order must be >= 0 (0 = constant conductivity).");
	}
	if (props.target_reflection <= 0.0 || props.target_reflection >= 1.0) {
		throw std::runtime_error("PML target_reflection must be in (0, 1).");
	}
	if (props.stretch_mode == PMLStretchMode::Box && props.radial_center.has_value()) {
		throw std::runtime_error(
			"PML radial_center is only valid when stretch_mode is \"radial\".");
	}

	return props;
}

} // namespace maxwell
