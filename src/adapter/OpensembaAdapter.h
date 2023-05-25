#pragma once

#include "components/Problem.h"
#include "solver/SolverOptions.h"

namespace maxwell {

class OpensembaAdapter {
	OpensembaAdapter(const std::string& filename);

	Problem readProblem() const;
	SolverOptions readSolverOptions() const;
};

}