#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include "components/probes.h"

using json = nlohmann::json;

namespace maxwell {

const FieldType assignFieldType(const std::string& field_type);

const Direction assignFieldSpatial(const std::string& direction);

Probes buildProbes(const json& case_data);

}