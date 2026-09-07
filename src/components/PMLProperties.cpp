#include "PMLProperties.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

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
	if (mat_json.contains("kappa_max") || mat_json.contains("alpha_max")) {
		throw std::runtime_error(
			"PML kappa_max/alpha_max are CFS-only and no longer accepted. "
			"See docs/pml/27-gedney-cfs-paused.md.");
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

	if (props.grading_order < 0) {
		throw std::runtime_error("PML grading_order must be >= 0 (0 = constant conductivity).");
	}
	if (props.target_reflection <= 0.0 || props.target_reflection >= 1.0) {
		throw std::runtime_error("PML target_reflection must be in (0, 1).");
	}

	return props;
}

} // namespace maxwell
