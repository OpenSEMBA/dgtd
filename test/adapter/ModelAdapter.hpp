#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include <TestUtils.h>

#include "components/model.h"

using json = nlohmann::json;

namespace maxwell {

using BoundaryPair = std::pair<AttributeToBoundary, AttributeToInteriorConditions>;

BdrCond assignBdrCond(const std::string& bdr_cond);

AttributeToMaterial assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh);
BoundaryPair        assembleAttributeToBoundary(const json& case_data, const mfem::Mesh& mesh);

std::string assembleMeshString(const std::string& filename);
mfem::Mesh assembleMesh(const std::string& mesh_string);

Model assembleModel(const json& case_data);

}