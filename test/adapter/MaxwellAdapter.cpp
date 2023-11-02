#pragma once
#include "MaxwellAdapter.hpp"

namespace maxwell {

Solver assembleCaseSolver(std::string case_name) 
{
	auto file_name{ maxwellCase(case_name) };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);
	Model model             { assembleModel        (case_data) };
	Probes probes           { assembleProbes       (case_data) };
	Sources sources         { assembleSources      (case_data) };
	SolverOptions solverOpts{ assembleSolverOptions(case_data) };

	return Solver(model, probes, sources, solverOpts);
}

}
