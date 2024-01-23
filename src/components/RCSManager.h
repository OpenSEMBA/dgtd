#pragma once

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/Probes.h>


namespace maxwell {

using namespace mfem;

using CompVec = std::vector<std::complex<double>>;
using FieldToRCSVector = std::map<std::string, CompVec>;
using FieldsToTime = std::pair<FieldToRCSVector, Time>;

class RCSManager {
public:

	RCSManager(const std::string& path, const std::string& probe_name, double f);

private:

	GridFunction getGridFunction(const std::string& path, const FieldType&, const Direction&);
	void performRCS2DCalculations(const std::string& path);

	Mesh m_;
	FieldToRCSVector fieldsRCS_;
	std::string basePath_;
	double frequency_{ 0.0 };

};

}
