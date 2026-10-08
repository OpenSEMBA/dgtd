#pragma once

#include <nlohmann/json.hpp>

#include "solver/Solver.h"

using json = nlohmann::json;

namespace maxwell::driver {
	json parseJSONfile(const std::string& case_name);

	mfem::Vector assemble3DVector(const json& input);

	maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest = true, bool restart = false);
	maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest, bool restart = false);

	std::string assembleMeshString(const std::string& filename);

	Probes buildProbes(const json& case_data);
	SolverOptions buildSolverOptions(const json& case_data);
	Sources buildSources(const json& case_data, const mfem::Mesh* mesh = nullptr);
	Model buildModel(
		const json& case_data,
		const std::string& case_path,
		const bool isTest,
		const int* partition_override = nullptr,
		int partition_count = 0);
}
