#pragma once
#include "ModelAdapter.hpp"

using FaceNo = int;
using isInterior = bool;

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
	else {
		throw std::exception(("The defined Boundary Type " + bdr_cond + " is incorrect.").c_str());
	}
}

std::string assembleMeshString(const std::string& filename)
{
	std::string folder_name{ filename };
	std::string s_msh  = ".msh";
	std::string s_mesh = ".mesh";

	std::string::size_type input_msh  = folder_name.find(s_msh);
	std::string::size_type input_mesh = folder_name.find(s_mesh);

	if (input_msh != std::string::npos)
		folder_name.erase(input_msh, s_msh.length());
	if (input_mesh != std::string::npos)
		folder_name.erase(input_mesh, s_mesh.length());

	return maxwellInputsFolder() + folder_name + "/" + filename;
}

AttributeToMaterial assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh)
{
	AttributeToMaterial res{};

	checkIfThrows(case_data.contains("model"), "JSON data does not include 'model'.");
	checkIfThrows(case_data["model"].contains("materials"), "JSON data does not include 'materials'.");

	auto show{ case_data["model"]["materials"] };

	for (auto m = 0; m < case_data["model"]["materials"].size(); m++) {
		for (auto t = 0; t < case_data["model"]["materials"][m]["tags"].size(); t++) {
			double eps{ 1.0 };
			double mu{ 1.0 };
			if (case_data["model"]["materials"][m].contains("relative_permittivity")) {
				eps = case_data["model"]["materials"][m]["relative_permittivity"];
			}
			if (case_data["model"]["materials"][m].contains("relative_permeability")) {
				mu = case_data["model"]["materials"][m]["relative_permeability"];
			}
			res.emplace(std::make_pair(case_data["model"]["materials"][m]["tags"][t], Material(eps, mu)));
		}
	}

	for (auto [att, v] : res) {
		checkIfThrows(
			mesh.attributes.Find(att) != -1, 
			std::string("There is no attribute") + std::to_string(att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}

	for (auto att{ 1 }; att < mesh.attributes.Size() + 1; att++) {
		checkIfThrows(
			!(res.find(att) == res.end()),
			std::string("There is no attribute") + std::to_string(att) +
			" defined in the JSON, but it is defined in the mesh."
		);
	}

	return res;
}

BoundaryPair assembleAttributeToBoundary(const json& case_data, const mfem::Mesh& mesh)
{
	checkIfThrows(case_data["model"].contains("boundaries"), 
		"JSON data does not include 'boundaries' in 'model'.");

	auto face2BdrEl{ mesh.GetFaceToBdrElMap() };

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {

		checkIfThrows(case_data["model"]["boundaries"][b].contains("tags"),
			"Boundary " + std::to_string(b) + " does not have defined 'tags'.");

		checkIfThrows(!case_data["model"]["boundaries"][b]["tags"].empty(),
			"Boundary " + std::to_string(b) + " 'tags' are empty.");

		checkIfThrows(case_data["model"]["boundaries"][b].contains("type"),
			"Boundary " + std::to_string(b) + " does not have a defined 'type'.");
	}

	std::map<Attribute, BdrCond>att2bdrCond;
	std::map<Attribute, isInterior> att2interior;
	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			att2bdrCond. emplace(case_data["model"]["boundaries"][b]["tags"][a], assignBdrCond(case_data["model"]["boundaries"][b]["type"]));
			att2interior.emplace(case_data["model"]["boundaries"][b]["tags"][a], false);
			for (auto f = 0; f < mesh.GetNumFaces(); f++) {
				auto bdrEl{ face2BdrEl[f] };
				if (bdrEl != -1) {
					if (mesh.GetBdrAttribute(bdrEl) == a && mesh.FaceIsInterior(f)) {
						att2interior[a] = true;
						break;
					}
				}
			}
		}
	}

	for (auto [att, v] : att2interior) {
		checkIfThrows(
			mesh.bdr_attributes.Find(att) != -1,
			std::string("There is no bdrattribute") + std::to_string(att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}

	for (auto att{ 1 }; att < mesh.bdr_attributes.Size() + 1; att++) {
		checkIfThrows(
			!(att2interior.find(att) == att2interior.end()),
			std::string("There is no bdrattribute") + std::to_string(att) +
			" defined in the JSON, but it is defined in the mesh."
		);
	}
	
	AttributeToBoundary att2bdr{};
	AttributeToInteriorConditions att2intCond{};
	for (auto [att, isInt] : att2interior) {
		switch (isInt) {
		case false:
			att2bdr.emplace(att, att2bdrCond[att]);
			break;
		case true:
			att2intCond.emplace(att, att2bdrCond[att]);
			break;
		}
	}
	
	return std::make_pair(att2bdr, att2intCond);
}

mfem::Mesh assembleMesh(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFile(mesh_string, 1, 0, true);
}

mfem::Mesh assembleMeshNoFix(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFileNoBdrFix(mesh_string, 1, 0, true);
}

Model assembleModel(const json& case_data)
{
	auto mesh{ assembleMesh(assembleMeshString(case_data["model"]["filename"]))};
		
	auto att_to_material{ assembleAttributeToMaterial(case_data, mesh) };
	auto att_to_bdr_info{ assembleAttributeToBoundary(case_data, mesh) };

	if (!att_to_bdr_info.second.empty()) {
		mesh = assembleMeshNoFix(assembleMeshString(case_data["model"]["filename"]));
	}

	return Model(mesh, att_to_material, att_to_bdr_info.first, att_to_bdr_info.second);
}
}