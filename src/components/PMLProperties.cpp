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
	if (token == "R" || token == "r") {
		throw std::runtime_error(
			"PML active_axes \"R\" is not supported. Use Cartesian \"X\"/\"Y\"/\"Z\" only.");
	}
	throw std::runtime_error(
		"PML active_axes entry must be \"X\", \"Y\", or \"Z\". Got: " + token);
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
			throw std::runtime_error(
				"PML stretch_mode \"radial\" is not supported. Use Cartesian box grading "
				"(omit stretch_mode, or set \"box\" / 0).");
		}
		throw std::runtime_error("PML stretch_mode integer must be 0 (box).");
	}
	if (v.is_string()) {
		std::string s = v.get<std::string>();
		std::transform(s.begin(), s.end(), s.begin(),
		               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		if (s == "box") {
			return PMLStretchMode::Box;
		}
		if (s == "radial") {
			throw std::runtime_error(
				"PML stretch_mode \"radial\" is not supported. Use Cartesian box grading "
				"(omit stretch_mode, or set \"box\" / 0).");
		}
		throw std::runtime_error(
			"PML stretch_mode must be \"box\" (or 0), or omitted. Got: " +
			v.get<std::string>());
	}
	throw std::runtime_error("PML stretch_mode must be a string or integer.");
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
	if (mat_json.contains("debye")) {
		throw std::runtime_error(
			"PML material must not define debye. Debye is a volumetric material outside the PML.");
	}
	if (mat_json.contains("matches_vacuum") && !mat_json["matches_vacuum"].get<bool>()) {
		throw std::runtime_error("Only matches_vacuum: true is supported for volumetric PML.");
	}
}

PMLProperties parsePMLMaterialBlock(const nlohmann::json& mat_json, int mesh_dim)
{
	validatePMLMaterialBlock(mat_json);

	PMLProperties props;
	props.matches_vacuum = mat_json.value("matches_vacuum", true);
	props.grading_order = mat_json.value("grading_order", 3);
	props.target_reflection = mat_json.value("target_reflection", 1e-6);
	props.stretch_mode = parseStretchMode(mat_json);
	props.kappa_max = mat_json.value("kappa_max", 1.0);
	props.alpha_max = mat_json.value("alpha_max", 0.0);
	props.active_axes = parseActiveAxes(mat_json, mesh_dim);

	if (mat_json.contains("radial_center")) {
		throw std::runtime_error(
			"PML radial_center is not supported. Cartesian box SC-PML grades from "
			"planar vacuum–PML interfaces.");
	}

	if (mat_json.contains("sigma_max")) {
		throw std::runtime_error(
			"PML sigma_max is not valid. Cartesian SC-PML uses target_reflection.");
	}
	if (props.alpha_max != 0.0) {
		throw std::runtime_error(
			"PML alpha_max > 0 is not supported for Cartesian SC-PML "
			"(CFS deferred). Use alpha_max: 0 or omit.");
	}
	if (props.target_reflection <= 0.0 || props.target_reflection >= 1.0) {
		throw std::runtime_error("PML target_reflection must be in (0, 1).");
	}

	if (props.kappa_max < 1.0) {
		throw std::runtime_error("PML kappa_max must be >= 1.");
	}
	if (props.grading_order < 0) {
		throw std::runtime_error("PML grading_order must be >= 0 (0 = constant conductivity).");
	}

	return props;
}

} // namespace maxwell
