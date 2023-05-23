#include <fstream>
#include <iostream>
#include <algorithm>

#include "ProblemDescription.h"

using namespace mfem;

namespace maxwell {

ProblemDescription ProblemDescription::readFromFile(const std::string& filename)
{
	std::ifstream f(filename);
	if (!f.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}
	
	ProblemDescription pD;
	// TODO stub

	return pD;
}

}
