#pragma once

#include "SourcesAdapter.hpp"

namespace maxwell {

using namespace fixtures::sources;

mfem::Vector assembleCenterVector(const json& source_center)
{
	mfem::Vector res(source_center.size());
	for (int i = 0; i < source_center.size(); i++) {
		res[i] = source_center[i];
	}
	return res;
}

mfem::Vector assemblePolarizationVector(const json& source_polarization)
{
	mfem::Vector res(3);
	for (int i = 0; i < source_polarization.size(); i++) {
		res[i] = source_polarization[i];
	}
	return res;
}

mfem::Vector assemblePropagationVector(const json& source_propagation)
{
	mfem::Vector res(3);
	for (int i = 0; i < source_propagation.size(); i++) {
		res[i] = source_propagation[i];
	}
	return res;
}

Sources assembleSources(const json& case_data)
{
	Sources res;
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] == "initial") {
			if (case_data["sources"][s]["magnitude"]["type"] == "gaussian") {
				return buildGaussianInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					case_data["sources"][s]["magnitude"]["spread"],
					assembleCenterVector(case_data["sources"][s]["center"]),
					assemblePolarizationVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["dimension"]
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "resonant") {
				return buildResonantModeInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemblePolarizationVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["magnitude"]["modes"]
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "totalField") {
			auto propagation(assemblePropagationVector(case_data));
			return buildGaussianPlanewave(
				case_data["sources"]["spread"],
				case_data["sources"]["delay"],
				assemblePolarizationVector(case_data["sources"][s]["polarization"]),
				assemblePropagationVector(case_data["sources"][s]["propagation"])
			);
		}
	}
}

}