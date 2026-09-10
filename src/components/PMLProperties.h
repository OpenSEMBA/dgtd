#pragma once

#include "Types.h"

#include <nlohmann/json.hpp>
#include <optional>
#include <set>
#include <vector>

namespace maxwell {

/// How σ depth is measured for classical ADE-PML grading.
/// Box: planar vacuum–PML interface per stretch axis (existing).
/// Radial: ρ = max(0, ‖x−c‖ − r_in); ADE stacks stay Cartesian (option A).
enum class PMLStretchMode {
	Box = 0,
	Radial = 1
};

struct PMLProperties {
	std::vector<Attribute> geom_tags;
	bool matches_vacuum = true;
	int grading_order = 3;
	double target_reflection = 1e-6;
	/// SC-PML stretch scale κ(ρ); default 1 recovers CuDG3D-equivalent ADE.
	double kappa_max = 1.0;
	/// CFS frequency shift (deferred; must remain 0 until pole is wired).
	double alpha_max = 0.0;
	std::set<Direction> active_axes;
	PMLStretchMode stretch_mode = PMLStretchMode::Box;
	/// Optional center for Radial mode. If unset, inferred from vacuum–PML interfaces.
	std::optional<std::array<double, 3>> radial_center;
};

std::set<Direction> parseActiveAxes(const nlohmann::json& mat_json, int mesh_dim);

void validatePMLMaterialBlock(const nlohmann::json& mat_json);

PMLProperties parsePMLMaterialBlock(const nlohmann::json& mat_json, int mesh_dim);

} // namespace maxwell
