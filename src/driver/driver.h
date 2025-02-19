#pragma once

#include <nlohmann/json.hpp>

#include "solver/Solver.h"

using json = nlohmann::json;

namespace maxwell::driver {
	json parseJSONfile(const std::string& case_name);

	maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest = true);
	maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest);

	std::string assembleMeshString(const std::string& filename);

	Probes buildProbes(const json& case_data);
	SolverOptions buildSolverOptions(const json& case_data);
	Sources buildSources(const json& case_data);
	Model buildModel(const json& case_data, const std::string& case_path, const bool isTest);
}