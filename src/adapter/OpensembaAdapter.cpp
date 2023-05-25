#include "OpensembaAdapter.h"

Problem Problem::readFromFile(const std::string& filename)
{
	std::ifstream f(filename);
	if (!f.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}

	Problem pD;
	// TODO stub

	return pD;
}