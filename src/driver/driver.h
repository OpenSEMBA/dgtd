#pragma once

#include <nlohmann/json.hpp>

#include "solver/Solver.h"
#include "components/Sources.h"
#include "math/Function.h"

using json = nlohmann::json;

namespace maxwell::driver {
	json parseJSONfile(const std::string& case_name);

	mfem::Vector assemble3DVector(const json& input);

	maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest = true);
	maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest);

	std::string assembleMeshString(const std::string& filename);

	Probes buildProbes(const json& case_data);
	SolverOptions buildSolverOptions(const json& case_data);
	Sources buildSources(const json& case_data);
	Model buildModel(const json& case_data, const std::string& case_path, const bool isTest);
}