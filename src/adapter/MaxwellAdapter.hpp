#include "solver/Solver.h"

#include "ModelAdapter.hpp"
#include "ProbesAdapter.hpp"
#include "SourcesAdapter.hpp"
#include "SolverOptsAdapter.hpp"

#include "../../test/TestUtils.h"

using json = nlohmann::json;

json parseJSONfile(const std::string& case_name);
maxwell::Solver buildSolver(const std::string& case_name);
maxwell::Solver buildSolver(const json& case_data);