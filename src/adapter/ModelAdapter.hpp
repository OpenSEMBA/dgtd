#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include "components/model.h"

#include "AdapterHelper.h"

using json = nlohmann::json;

using FaceNo = int;
using isInterior = bool;

namespace maxwell {

using BoundaryPair = std::pair<GeomTagToBoundary, GeomTagToInteriorConditions>;

BdrCond assignBdrCond(const std::string& bdr_cond);

GeomTagToMaterial assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh);
BoundaryPair      assembleAttributeToBoundary(const json& case_data, const mfem::Mesh& mesh);

std::string assembleMeshString(const std::string& filename);
mfem::Mesh assembleMesh(const std::string& mesh_string);

Model buildModel(const json& case_data);

}