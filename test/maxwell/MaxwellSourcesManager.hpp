#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include "components/sources.h"
#include "SourceFixtures.h"
#include "MaxwellProbesManager.hpp"

using json = nlohmann::json;

namespace maxwell {

using namespace fixtures::sources;

mfem::Vector assembleCenterVector(const json& case_data);

mfem::Vector assemblePolarizationVector(const json& case_data);

Source::CartesianAngles assembleRotationVector(const json& case_data);

mfem::Vector assemblePropagationVector(const json& case_data);

Sources assembleSources(const json& case_data);

}