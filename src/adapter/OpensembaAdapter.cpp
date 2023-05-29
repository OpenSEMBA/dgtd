#include "OpensembaAdapter.h"
#include "core/ProblemDescription.h"

namespace maxwell {

	OpensembaAdapter::OpensembaAdapter(const std::string& fn) :
	filename_{fn}
{
	std::ifstream f(filename_);
	if (!f.is_open()) {
		throw std::runtime_error("Could not open file: " + filename_);
	}

	SEMBA::UnstructuredProblemDescription pd;
}

Problem OpensembaAdapter::readProblem() const
{
	
	Problem r;
	// TODO stub

	return r;
}

SolverOptions OpensembaAdapter::readSolverOptions() const
{

	SolverOptions r;
	// TODO stub

	return r;
}

}