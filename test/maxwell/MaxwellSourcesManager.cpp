#pragma once

#include "MaxwellSourcesManager.hpp"

namespace maxwell {

using namespace fixtures::sources;

mfem::Vector assembleCenterVector(const json& case_data)
{
	mfem::Vector res(case_data["dimension"]);
	for (int i = 0; i < case_data["dimension"]; i++) {
		res[i] = case_data["sources"]["center"][i];
	}
	return res;
}

mfem::Vector assemblePolarizationVector(const json& case_data)
{
	mfem::Vector res(case_data["sources"]["polarization"]);
	for (int i = 0; i < case_data["sources"]["polarization"].size(); i++) {
		res[i] = case_data["sources"]["polarization"][i];
	}
	return res;
}

Source::CartesianAngles assembleRotationVector(const json& case_data)
{
	Source::CartesianAngles res(3);
	res[0] = 0.0;
	res[1] = 0.0;
	res[2] = 0.0;
	if (case_data["rotation_angles"]) {
		for (int i = 0; i < case_data["sources"]["rotation_angles"].size(); i++) {
			res[i] = case_data["sources"]["rotation_angles"][i];
		}
	}
	return res;
}

mfem::Vector assemblePropagationVector(const json& case_data)
{
	mfem::Vector res(case_data["sources"]["propagation"]);
	for (int i = 0; i < case_data["sources"]["propagation"].size(); i++) {
		res[i] = case_data["sources"]["propagation"][i];
	}
	return res;
}

Sources assembleSources(const json& case_data)
{
	
	auto field{ assignFieldType(case_data["sources"]["field_type"]) };
	auto direction{ assignFieldSpatial(case_data["sources"]["field_spatial"])};
	auto center(assembleCenterVector(case_data));
	auto polarization(assemblePolarizationVector(case_data));
	auto rotation_angles(assembleRotationVector(case_data));

	Sources res;
	if (case_data["sources"]["type"] == "Initial") {
		return buildGaussianInitialField(
			field, 
			case_data["sources"]["spread"], 
			center, 
			polarization, 
			case_data["sources"]["source_dimension"], 
			rotation_angles
		);
	}
	else if (case_data["sources"]["type"] == "Resonant") {
		return buildResonantModeInitialField(
			field,
			polarization,
			case_data["sources"]["modes"]
		);
	}
	else if (case_data["sources"]["type"] == "TDPlanewave") {
		auto propagation(assemblePropagationVector(case_data));
		return buildGaussianPlanewave(
			case_data["sources"]["spread"],
			case_data["sources"]["delay"],
			polarization,
			propagation
		);
	}
	else if (case_data["sources"]["type"] == "InitPlanewave") {
		auto propagation(assemblePropagationVector(case_data));
		if (case_data["sources"]["function_type"] == "Gaussian") {
			auto function{ 
				Gaussian(
					case_data["sources"]["spread"], 
					mfem::Vector(0.0),
					case_data["sources"]["sources_dimension"]
				)
			};
			return buildPlanewaveInitialField(
				function, 
				center, 
				polarization, 
				propagation, 
				rotation_angles
			);
		}
	}


}

}