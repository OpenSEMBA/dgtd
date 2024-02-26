#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include <components/Sources.h>
#include <math/SourceFixtures.h>
#include <adapter/ProbesAdapter.hpp>

using json = nlohmann::json;

namespace maxwell {

using namespace fixtures::sources;

mfem::Vector assembleCenterVector(const json& case_data);

mfem::Vector assemble3DVector(const json& case_data);

Sources buildSources(const json& case_data);

}