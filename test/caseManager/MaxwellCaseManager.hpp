#include "solver/Solver.h"

#include "MaxwellModelManager.hpp"
#include "MaxwellProbesManager.hpp"
#include "MaxwellSourcesManager.hpp"
#include "MaxwellSolverOptsManager.hpp"

using json = nlohmann::json;

maxwell::Solver assembleCaseSolver(std::string case_name);