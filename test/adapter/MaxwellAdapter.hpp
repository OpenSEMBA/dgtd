#include "solver/Solver.h"

#include "ModelAdapter.hpp"
#include "ProbesAdapter.hpp"
#include "SourcesAdapter.hpp"
#include "SolverOptsAdapter.hpp"

using json = nlohmann::json;

maxwell::Solver assembleCaseSolver(std::string case_name);