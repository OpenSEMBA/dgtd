#pragma once

#include "Types.h"

#include <nlohmann/json.hpp>
#include <set>
#include <vector>

namespace maxwell {

struct PMLProperties {
	std::vector<Attribute> geom_tags;
	bool matches_vacuum = true;
	int grading_order = 3;
	double target_reflection = 1e-6;
	double kappa_max = 1.0;
	double alpha_max = 0.0;
	std::set<Direction> active_axes;
};

std::set<Direction> parseActiveAxes(const nlohmann::json& mat_json, int mesh_dim);

PMLProperties parsePMLMaterialBlock(const nlohmann::json& mat_json, int mesh_dim);

void validatePMLMaterialBlock(const nlohmann::json& mat_json);

} // namespace maxwell
