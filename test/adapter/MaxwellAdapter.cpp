#pragma once
#include "MaxwellAdapter.hpp"

json parseJSONfile(const std::string& case_name)
{
	auto file_name{ maxwellCase(case_name) };
	std::ifstream test_file(file_name);
	return json::parse(test_file);
}

maxwell::Solver buildSolver(const std::string& case_name)
{
	auto case_data{ parseJSONfile(case_name) };

	return buildSolver(case_data);
}

void postProcessInformation(const json& case_data, maxwell::Model& model) 
{
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
	mfem::Array<int> tfsf_tags;
		if (case_data["sources"][s]["type"] == "totalField") {
			for (auto t{ 0 }; t < case_data["sources"][s]["tags"].size(); t++) {
				tfsf_tags.Append(case_data["sources"][s]["tags"][t]);
			}
		auto marker{ model.getMarker(maxwell::BdrCond::TotalFieldIn, true) };
		marker.SetSize(model.getConstMesh().bdr_attributes.Max());
		marker = 0;
		for (auto t : tfsf_tags) {
			marker[t - 1] = 1;
		}
		model.getTotalFieldScatteredFieldToMarker().insert(std::make_pair(maxwell::BdrCond::TotalFieldIn, marker));
		}
	}
}

maxwell::Solver buildSolver(const json& case_data)
{
	maxwell::Model model{ maxwell::buildModel(case_data) };
	maxwell::Probes probes{ maxwell::buildProbes(case_data) };
	maxwell::Sources sources{ maxwell::buildSources(case_data) };
	maxwell::SolverOptions solverOpts{ maxwell::buildSolverOptions(case_data) };

	postProcessInformation(case_data, model);

	return maxwell::Solver(model, probes, sources, solverOpts);
}
