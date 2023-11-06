#pragma once
#include "MaxwellAdapter.hpp"

maxwell::Solver assembleCaseSolver(std::string case_name)
{
	auto file_name{ maxwellCase(case_name) };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);
	maxwell::Model model{ maxwell::assembleModel(case_data) };
	maxwell::Probes probes{ maxwell::assembleProbes(case_data) };
	maxwell::Sources sources{ maxwell::assembleSources(case_data) };
	maxwell::SolverOptions solverOpts{ maxwell::assembleSolverOptions(case_data) };

	return maxwell::Solver(model, probes, sources, solverOpts);
}
