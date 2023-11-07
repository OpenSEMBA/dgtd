#pragma once
#include "MaxwellAdapter.hpp"

json parseJSONfile(const std::string& case_name)
{
	auto file_name{ maxwellCase(case_name) };
	std::ifstream test_file(file_name);
	return json::parse(test_file);
}

maxwell::Solver buildSolver(const std::string& case_name)
{
	auto case_data{ parseJSONfile(case_name) };

	return buildSolver(case_data);
}

maxwell::Solver buildSolver(const json& case_data)
{
	maxwell::Model model{ maxwell::buildModel(case_data) };
	maxwell::Probes probes{ maxwell::buildProbes(case_data) };
	maxwell::Sources sources{ maxwell::buildSources(case_data) };
	maxwell::SolverOptions solverOpts{ maxwell::buildSolverOptions(case_data) };

	return maxwell::Solver(model, probes, sources, solverOpts);
}
