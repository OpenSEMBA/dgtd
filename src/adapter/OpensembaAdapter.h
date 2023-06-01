#pragma once

#include <nlohmann/json.hpp>

#include "components/Problem.h"
#include "solver/SolverOptions.h"

namespace maxwell {

class OpensembaAdapter {
public:
	OpensembaAdapter(const std::string& filename);

	Problem readProblem() const;
	SolverOptions readSolverOptions() const;

private:
	std::string filename_;
	nlohmann::json json_;
};

}