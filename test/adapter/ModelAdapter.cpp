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

GeomTagToMaterial assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh)
{
	GeomTagToMaterial res{};

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

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {

		checkIfThrows(case_data["model"]["boundaries"][b].contains("tags"),
			"Boundary " + std::to_string(b) + " does not have defined 'tags'.");

		checkIfThrows(!case_data["model"]["boundaries"][b]["tags"].empty(),
			"Boundary " + std::to_string(b) + " 'tags' are empty.");

		checkIfThrows(case_data["model"]["boundaries"][b].contains("class"),
			"Boundary " + std::to_string(b) + " does not have a defined 'class'.");
	}

	auto face2BdrEl{ mesh.GetFaceToBdrElMap() };

	std::map<GeomTag, BdrCond> geomTag2bdrCond;
	std::map<GeomTag, isInterior> geomTag2interior;
	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			geomTag2bdrCond. emplace(case_data["model"]["boundaries"][b]["tags"][a], assignBdrCond(case_data["model"]["boundaries"][b]["class"]));
			geomTag2interior.emplace(case_data["model"]["boundaries"][b]["tags"][a], false);
			for (auto f = 0; f < mesh.GetNumFaces(); f++) {
				if (face2BdrEl[f] != -1) {
					if (mesh.GetBdrAttribute(face2BdrEl[f]) == a && mesh.FaceIsInterior(f)) {
						geomTag2interior[a] = true;
						break;
					}
				}
			}
		}
	}

	for (auto [geomTag, v] : geomTag2interior) {
		checkIfThrows(
			mesh.bdr_attributes.Find(geomTag) != -1,
			std::string("There is no boundary geometrical tag ") + std::to_string(geomTag) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}

	for (auto geomTag{ 1 }; geomTag < mesh.bdr_attributes.Size() + 1; geomTag++) {
		checkIfThrows(
			!(geomTag2interior.find(geomTag) == geomTag2interior.end()),
			std::string("There is no  boundary geometrical tag ") + std::to_string(geomTag) +
			" defined in the JSON, but it is defined in the mesh."
		);
	}
	
	GeomTagToBoundary geomTag2bdr{};
	GeomTagToInteriorConditions geomTag2intCond{};
	for (auto [att, isInt] : geomTag2interior) {
		switch (isInt) {
		case false:
			geomTag2bdr.emplace(att, geomTag2bdrCond[att]);
			break;
		case true:
			geomTag2intCond.emplace(att, geomTag2bdrCond[att]);
			break;
		}
	}
	
	return std::make_pair(geomTag2bdr, geomTag2intCond);
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