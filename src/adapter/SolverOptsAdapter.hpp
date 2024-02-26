#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include <solver/SolverOptions.h>

using json = nlohmann::json;

namespace maxwell {

SolverOptions buildSolverOptions(const json& case_data);

}
