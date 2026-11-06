#pragma once

#include <nlohmann/json.hpp>

#include "solver/Solver.h"

using json = nlohmann::json;

namespace maxwell::driver {
	json parseJSONfile(const std::string& case_name);

	mfem::Vector assemble3DVector(const json& input);

	/// Gaussian σ from `f_1e` (Hz, 1/e incident power) or `spread` (light-metres).
	/// If both are set, uses `f_1e` and warns. Throws if neither is set.
	double assembleGaussianSpread(const json& obj, const char* context);

	enum class LaunchResume {
		None,
		Continue,
		Restart
	};

	maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest = true, LaunchResume resume = LaunchResume::None);
	maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest, LaunchResume resume = LaunchResume::None);

	std::string assembleMeshString(const std::string& filename);

	Probes buildProbes(const json& case_data);
	SolverOptions buildSolverOptions(const json& case_data);
	Sources buildSources(const json& case_data, const mfem::Mesh* mesh = nullptr, const std::string& case_path = "");
	Model buildModel(
		const json& case_data,
		const std::string& case_path,
		const bool isTest,
		const int* partition_override = nullptr,
		int partition_count = 0);
}
