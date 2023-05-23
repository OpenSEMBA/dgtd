#pragma once

#include "Types.h"
#include "Fields.h"
#include "ProbesManager.h"
#include "SourcesManager.h"
#include "SolverOptions.h"
#include "evolution/Evolution3D.h"
#include "evolution/Evolution3D_Spectral.h"

namespace maxwell {

struct ProblemDescription {
    Model model;
    Probes probes;
    Sources sources;

    static ProblemDescription readFromFile(const std::string& filename);
};

}