#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include <TestUtils.h>

#include "components/model.h"

using json = nlohmann::json;

namespace maxwell {

BdrCond assignBdrCond(const std::string& bdr_cond);
BdrCond assignInteriorCond(const std::string& bdr_cond);

mfem::Element::Type assignElementType(const std::string& element_type);
mfem::Element::Type assignElementType(const std::string& element_type, const int dimension);

mfem::Vector assembleNumberOfElements(const json& case_data);

mfem::Element::Type assembleElementType(const json& case_data);

mfem::Vector assembleMeshLength(const json& case_data);

void verifyDimensionAndElementCompatibility(const int dimension, const mfem::Element::Type& element_type);

std::string assembleMeshString(const std::string& mesh_name, const std::string& mesh_format);

AttributeToMaterial assembleAttributeToMaterial(const json& case_data);

std::pair<AttributeToBoundary, bool> assembleAttributeToBoundary(const json& case_data);

AttributeToInteriorConditions assembleAttributeToInteriorConditions(const json& case_data);

mfem::Mesh assembleAutoMesh(
	const mfem::Vector& num_of_elem,
	const mfem::Element::Type,
	const mfem::Vector& mesh_length);

mfem::Mesh assembleManualMesh(const std::string& mesh_string, bool totalfieldscatteredfield_flag);

Model assembleModel(const json& case_data);
}