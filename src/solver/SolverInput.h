#pragma once

#include "components/Problem.h"
#include "SolverOptions.h"

namespace maxwell {

struct SolverInput {
    Problem problem;
    SolverOptions options;
};

}
