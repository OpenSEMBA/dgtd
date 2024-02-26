#pragma once
#include "ProbesAdapter.hpp"

namespace maxwell {

const FieldType assignFieldType(const std::string& field_type)
{
	if (field_type == "E") {
		return FieldType::E;
	}
	else if (field_type == "H") {
		return FieldType::H;
	}
	else {
		throw std::exception("Wrong Field Type in Point Probe assignation.");
	}
}

const Direction assignFieldSpatial(const std::string& direction)
{
	if (direction == "X") {
		return X;
	}
	else if (direction == "Y") {
		return Y;
	}
	else if (direction == "Z") {
		return Z;
	}
	else {
		throw std::exception("Wrong Field Polarization in Point Probe assignation.");
	}
}

std::vector<double> assembleVector(const json& input)
{
	std::vector<double> res(input.size());
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
}

Probes buildProbes(const json& case_data) 
{

	Probes probes;
	
	if (case_data["probes"].contains("exporter")) {
		ExporterProbe exporter_probe;
		exporter_probe.name = case_data["model"]["filename"];
		if (case_data["probes"]["exporter"].contains("steps")) {
			exporter_probe.visSteps = case_data["probes"]["exporter"]["steps"];
		}
		probes.exporterProbes.push_back(exporter_probe);
	}

	if (case_data["probes"].contains("field")) {
		for (int p = 0; p < case_data["probes"]["field"].size(); p++) {
			FieldProbe field_probe(
				assembleVector(case_data["probes"]["field"][p]["position"])
			);
			probes.fieldProbes.push_back(field_probe);
		}
	}

	if (case_data["probes"].contains("neartofarfield")) {
		for (int p{ 0 }; p < case_data["probes"]["neartofarfield"].size(); p++) {
			NearToFarFieldProbe probe;
			if (case_data["probes"]["neartofarfield"][p].contains("name")) {
				probe.name = case_data["probes"]["neartofarfield"][p]["name"];
			}
			if (case_data["probes"]["neartofarfield"][p].contains("steps")) {
				probe.steps = case_data["probes"]["neartofarfield"][p]["steps"];
			}
			if (case_data["probes"]["neartofarfield"][p].contains("tags")) {
				std::vector<int> tags;
				for (int t{ 0 }; t < case_data["probes"]["neartofarfield"][p]["tags"].size(); t++) {
					tags.push_back(case_data["probes"]["neartofarfield"][p]["tags"][t]);
				}
				probe.tags = tags;
			}
			else {
				throw std::exception("Tags have not been defined in neartofarfield probe.");
			}
			probes.nearToFarFieldProbes.push_back(probe);
		}
	}

	return probes;
}

}