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

Probes assembleProbes(const json& case_data) 
{

	Probes probes;
	
	if (case_data["exporter_probes"]) {
		ExporterProbe exporter_probe;
		exporter_probe.name = case_data["name"];
		exporter_probe.visSteps = case_data["exporter_probe"]["steps"];
		probes.exporterProbes.push_back(exporter_probe);
	}

	if (case_data["point_probes"]) {
		for (int p = 0; p < case_data["point_probes"].size(); p++) {
			auto field{ assignFieldType(case_data["point_probes"][p][0]) };
			auto direction{ assignFieldSpatial(case_data["point_probes"][p][1]) };
			std::vector<double> point(case_data["point_probes"][p][2].size());
			for (int i = 0; i < point.size(); i++) {
				point[i] = case_data["point_probes"][p][2][i];
			}
			PointProbe point_probe(field, direction, point);
			probes.pointProbes.push_back(point_probe);
		}
	}
	return probes;
}

}