#include <fstream>
#include <nlohmann/json.hpp>

#include <TestUtils.h>

#include "components/model.h"

using json = nlohmann::json;

namespace maxwell {

BdrCond assignBdrCond(const std::string& bdr_cond)
{
	if (bdr_cond == "PEC") {
		return BdrCond::PEC;
	}
	else if (bdr_cond == "PMC") {
		return BdrCond::PMC;
	}
	else if (bdr_cond == "SMA") {
		return BdrCond::SMA;
	}
	else if (bdr_cond == "NONE") {
		return BdrCond::NONE;
	}
	else if (bdr_cond == "TotalFieldIn") {
		return BdrCond::TotalFieldIn;
	}
}

BdrCond assignInteriorCond(const std::string& bdr_cond)
{
	return assignBdrCond(bdr_cond);
}

mfem::Element::Type assignElementType(const std::string& element_type)
{
	if (element_type == "Triangle") {
		return mfem::Element::Type::TRIANGLE;
	}
	else if (element_type == "Quadrilateral") {
		return mfem::Element::Type::QUADRILATERAL;
	}
	else if (element_type == "Tetrahedron") {
		return mfem::Element::Type::TETRAHEDRON;
	}
	else if (element_type == "Pyramid") {
		return mfem::Element::Type::PYRAMID;
	}
	else if (element_type == "Wedge") {
		return mfem::Element::Type::WEDGE;
	}
	else if (element_type == "Hexahedron") {
		return mfem::Element::Type::HEXAHEDRON;
	}
}

mfem::Element::Type assignElementType(const std::string& element_type, const int dimension)
{
	assert(element_type == "Undefined");
	switch (dimension) {
	case 1:
		return mfem::Element::Type::SEGMENT;
	case 2:
		return mfem::Element::Type::TRIANGLE;
	case 3:
		return mfem::Element::Type::TETRAHEDRON;
	default:
		throw std::exception("Incorrect dimension assigning element type for case.");
	}
}

mfem::Vector assembleNumberOfElements(const json& case_data)
{
	mfem::Vector res(case_data["dimension"]);
	if (case_data.contains("number_of_elements")) {
		for (int i = 0; i < case_data["dimension"]; i++) {
			res[i] = case_data["number_of_elements"][i];
		}
	}
	return res;
}

mfem::Element::Type assembleElementType(const json& case_data)
{
	if (case_data.contains("element_type")) {
		return assignElementType(case_data["element_type"]);
	}
	else {
		return assignElementType("Undefined", case_data["dimension"]);
	}
}

mfem::Vector assembleMeshLength(const json& case_data) {
	mfem::Vector res(case_data["dimension"]);
	res = 1.0;
	if (case_data.contains("mesh_length")) {
		for (int i = 0; i < case_data["dimension"]; i++) {
			res[i] = case_data["mesh_length"][i];
		}
	}
	return res;
}

void verifyDimensionAndElementCompatibility(const int dimension, const mfem::Element::Type& element_type)
{
	switch (dimension) {
	case 2:
		if (element_type == mfem::Element::Type::TRIANGLE
			|| element_type == mfem::Element::Type::QUADRILATERAL) {
			break;
		}
		else {
			throw std::exception("Incompatible dimension and element type for case.");
		}
	case 3:
		if (element_type == mfem::Element::Type::TETRAHEDRON
			|| element_type == mfem::Element::Type::PYRAMID
			|| element_type == mfem::Element::Type::WEDGE
			|| element_type == mfem::Element::Type::HEXAHEDRON) {
			break;
		}
		else {
			throw std::exception("Incompatible dimension and element type in case.");
		}
	default:
		break;

	}
}

std::string assembleMeshString(const std::string& mesh_name, const std::string& mesh_format)
{
	std::string case_name{ mesh_name };
	std::string m_format;
	if (mesh_format == "gmsh") {
		m_format = ".msh";
	}
	else if (mesh_format == "mfem") {
		m_format = ".mesh";
	}

	return maxwellInputsFolder() + case_name + "/" + mesh_name + m_format;
}

AttributeToMaterial assembleAttributeToMaterial(const json& case_data)
{
	AttributeToMaterial res{};
	if (case_data.contains("attribute_to_material")) {
		for (int i = 0; i < case_data["attribute_to_material"].size(); i++) {
			Material mat(
				case_data["attribute_to_material"][i][1][0], //Epsilon
				case_data["attribute_to_material"][i][1][1]  //Mu
			);
			res.emplace(std::make_pair(case_data["attribute_to_material"][i][0], mat));
		}
	}
	return res;
}

std::pair<AttributeToBoundary, bool> assembleAttributeToBoundary(const json& case_data)
{
	AttributeToBoundary att_to_boundary{};
	bool totalfield_flag = false;
	if (case_data.contains("attribute_to_boundary")) {
		for (int i = 0; i < case_data["attribute_to_boundary"].size(); i++) {
			std::pair<Attribute, BdrCond> pair;
			pair.first = case_data["attribute_to_boundary"][i][0];
			pair.second = assignBdrCond(case_data["attribute_to_boundary"][i][1]);
			att_to_boundary.emplace(pair);
			if (pair.first == 301) {
				totalfield_flag = true;
			}
		}
	}
	return std::make_pair(att_to_boundary, totalfield_flag);
}

AttributeToInteriorConditions assembleAttributeToInteriorConditions(const json& case_data)
{
	AttributeToInteriorConditions res{};
	if (case_data.contains("attribute_to_interior_conditions")) {
		for (int i = 0; i < case_data["attribute_to_interior_conditions"].size(); i++) {
			std::pair<Attribute, BdrCond> pair;
			pair.first = case_data["attribute_to_interior_conditions"][i][0];
			pair.second = assignBdrCond(case_data["attribute_to_interior_conditions"][i][1]);
			res.emplace(pair);
		}
	}
	return res;
}

mfem::Mesh assembleAutoMesh(
	const mfem::Vector& num_of_elem,
	const mfem::Element::Type element_type,
	const mfem::Vector& mesh_length)
{
	mfem::Mesh mesh;
	switch (num_of_elem.Size()) {
	case 1:
		mesh = mfem::Mesh::MakeCartesian1D(
			num_of_elem[0],
			mesh_length[0]);
		break;
	case 2:
		mesh = mfem::Mesh::MakeCartesian2D(
			num_of_elem[0], num_of_elem[1],
			element_type, true,
			mesh_length[0], mesh_length[1]);
		break;
	case 3:
		mesh = mfem::Mesh::MakeCartesian3D(
			num_of_elem[0], num_of_elem[1], num_of_elem[2],
			element_type,
			mesh_length[0], mesh_length[1], mesh_length[2]);
		break;
	default:
		throw std::exception("Incorrect dimension in assembleModel for case.");
	}
	return mesh;
}

mfem::Mesh assembleManualMesh(const std::string& mesh_string, bool totalfieldscatteredfield_flag)
{
	switch (totalfieldscatteredfield_flag) {
	case true:
		return mfem::Mesh::LoadFromFileNoBdrFix(mesh_string, 1, 0, true);
	case false:
		return mfem::Mesh::LoadFromFile(mesh_string, 1, 0, true);
	}
}

Model assembleModel(const json& case_data)
{
	if (case_data.contains("auto_model")) {

		auto num_of_elem { assembleNumberOfElements(case_data) };
		auto element_type{ assembleElementType     (case_data) };
		auto mesh_length { assembleMeshLength      (case_data) };

		verifyDimensionAndElementCompatibility(case_data["dimension"], element_type);

		auto mesh{ assembleAutoMesh(num_of_elem, element_type, mesh_length) };

		AttributeToMaterial att_to_material{ assembleAttributeToMaterial(case_data) };
		auto bdr_data{ assembleAttributeToBoundary(case_data) };
		AttributeToBoundary att_to_boundary{ bdr_data.first };
		AttributeToInteriorConditions att_to_int_conds{ assembleAttributeToInteriorConditions(case_data) };

		return Model(mesh, att_to_material, att_to_boundary, att_to_int_conds);
	}
	else if (case_data.contains("manual_model")) {

		AttributeToMaterial att_to_material{ assembleAttributeToMaterial(case_data) };

		auto bdr_data{ assembleAttributeToBoundary(case_data) };
		AttributeToBoundary att_to_boundary{ bdr_data.first };
		bool tfsf_flag{ bdr_data.second };
		AttributeToInteriorConditions att_to_int_conds{ assembleAttributeToInteriorConditions(case_data) };

		auto mesh{ assembleManualMesh(
			assembleMeshString(
				case_data["mesh"],
				case_data["mesh_format"]),
			tfsf_flag)
		};

		return Model(mesh, att_to_material, att_to_boundary, att_to_int_conds);
	}
	else {
		throw std::exception("Json File does not contain a valid model input.");
	}
}
}