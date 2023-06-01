#pragma once

#include <nlohmann/json.hpp>

#include "components/Problem.h"
#include "solver/SolverOptions.h"

namespace maxwell {

class OpensembaAdapter {
public:
	OpensembaAdapter(const std::string& filename);

	Problem readProblem() const;
	SolverOptions readSolverOptions() const;

private:
	using json = nlohmann::json;

	std::string filename_;
	json json_;

	Model readModel(const json& j) const;
	mfem::Mesh readMesh(const json& j) const;
};

}