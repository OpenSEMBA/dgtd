#include <fstream>
#include <nlohmann/json.hpp>

#include <TestUtils.h>

#include "components/sources.h"
#include "MaxwellProbesManager.cpp"

using json = nlohmann::json;

namespace maxwell {


Sources assembleSources(const json& case_data)
{
	Sources res;
	auto field{ assignFieldType(case_data["sources"]["field_type"]) };
	auto direction{ assignFieldSpatial(case_data["sources"]["field_spatial"])};
	mfem::Vector center(case_data["dimension"]);
	for (int i = 0; i < case_data["dimension"]; i++) {
		center[i] = case_data["sources"]["center"][i];
	}
	mfem::Vector polarization(case_data["sources"]["polarization"]);
	for (int i = 0; i < case_data["sources"]["polarization"].size(); i++) {
		polarization[i] = case_data["sources"]["polarization"][i];
	}
	const int source_dimension{ case_data["source_dimension"] };
	mfem::Vector rotation_angles(3);
	rotation_angles = 0.0;
	if (case_data["rotation_angles"]) {
		for (int i = 0; i < case_data["sources"]["rotation_angles"].size(); i++) {
			rotation_angles = case_data["sources"]["rotation_angles"][i];
		}
	}

	return res;
}

}