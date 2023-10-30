#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include <solver/SolverOptions.h>

using json = nlohmann::json;

namespace maxwell {

SolverOptions assembleSolverOptions(const json& case_data) 
{
	SolverOptions res;
	if (case_data["solver_options"]["solver_type"] == "Centered") {
		res.setCentered();
	}
	else if (case_data["solver_options"]["solver_type"] == "Upwind") {}
	if (case_data["solver_options"]["time_step"]) {
		res.setTimeStep(case_data["solver_options"]["time_step"]);
	}
	if (case_data["solver_options"]["final_time"]) {
		res.setFinalTime(case_data["solver_options"]["final_time"]);
	}
	if (case_data["solver_options"]["cfl"]) {
		res.setCFL(case_data["solver_options"]["cfl"]);
	}
	if (case_data["solver_options"]["order"]) {
		res.setOrder(case_data["solver_options"]["order"]);
	}
	if (case_data["solver_options"]["spectral"]) {
		res.setSpectralEO(case_data["solver_options"]["spectral"]);
	}
	return res;
}

}
