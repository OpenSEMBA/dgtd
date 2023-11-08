#pragma once

#include "SourcesAdapter.hpp"

namespace maxwell {

using namespace fixtures::sources;

mfem::Vector assembleCenterVector(const json& source_center)
{
	mfem::Vector res(int(source_center.size()));
	for (int i = 0; i < source_center.size(); i++) {
		res[i] = source_center[i];
	}
	return res;
}

mfem::Vector assemble3DVector(const json& input)
{
	mfem::Vector res(3);
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
}

Sources buildSources(const json& case_data)
{
	Sources res;
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] == "initial") {
			if (case_data["sources"][s]["magnitude"]["type"] == "gaussian") {
				return buildGaussianInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					case_data["sources"][s]["magnitude"]["spread"],
					assembleCenterVector(case_data["sources"][s]["center"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["dimension"]
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "resonant") {
				return buildResonantModeInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["magnitude"]["modes"]
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "totalField") {
			return buildGaussianPlanewave(
				case_data["sources"][s]["magnitude"]["spread"],
				case_data["sources"][s]["magnitude"]["delay"],
				assemble3DVector(case_data["sources"][s]["polarization"]),
				assemble3DVector(case_data["sources"][s]["propagation"])
			);
		}
		else {
			throw std::exception("Unknown source type in Json.");
		}
	}
}

}