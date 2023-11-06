#pragma once
#include "SolverOptsAdapter.hpp"

namespace maxwell {

SolverOptions assembleSolverOptions(const json& case_data) 
{
	SolverOptions res{};

	if (case_data.contains("solver_options")) {

		if (case_data["solver_options"].contains("solver_type")) {
			if (case_data["solver_options"]["solver_type"] == "centered")
				res.setCentered();
			else if (case_data["solver_options"]["solver_type"] == "upwind") {}
		}

		if (case_data["solver_options"].contains("time_step")) {
			res.setTimeStep(case_data["solver_options"]["time_step"]);
		}

		if (case_data["solver_options"].contains("final_time")) {
			res.setFinalTime(case_data["solver_options"]["final_time"]);
		}

		if (case_data["solver_options"].contains("cfl")) {
			res.setCFL(case_data["solver_options"]["cfl"]);
		}

		if (case_data["solver_options"].contains("order")) {
			res.setOrder(case_data["solver_options"]["order"]);
		}

		if (case_data["solver_options"].contains("spectral")) {
			res.setSpectralEO(case_data["solver_options"]["spectral"]);
		}

	}
	return res;
}

}
