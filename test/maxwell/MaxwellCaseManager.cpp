#include <fstream>
#include <nlohmann/json.hpp>

#include <TestUtils.h>

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"

#include <MaxwellModelManager.cpp>
#include <MaxwellProbesManager.cpp>
#include <MaxwellSourcesManager.cpp>

using json = nlohmann::json;

namespace maxwell {

Solver assembleCaseSolver(std::string case_name) 
{
	auto file_name{ maxwellCase("1D_PEC_Centered") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);
	Model model{ assembleModel(case_data) };
	Probes probes{ assembleProbes(case_data) };
	Sources sources{ assembleSources(case_data) };
}

}
