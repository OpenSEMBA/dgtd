#pragma once

#include "Model.h"
#include "Probes.h"
#include "Sources.h"

namespace maxwell {

struct ProblemDescription {
    Model model;
    Probes probes;
    Sources sources;
};

}