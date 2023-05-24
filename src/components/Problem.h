#pragma once

#include "Model.h"
#include "Probes.h"
#include "Sources.h"

namespace maxwell {

struct Problem {
    Model model;
    Probes probes;
    Sources sources;
};

}