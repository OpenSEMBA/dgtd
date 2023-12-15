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


	RCSManager(const std::string& path, const NearToFarFieldProbe&);

	void update(FieldsToTime&);

private:

	void initFieldsRCS(const std::string&);
	void calculateRCS(CompVec&, const GridFunction&, const Time);
	GridFunction getGridFunction(const std::string& path, const FieldType&, const Direction&);

	std::unique_ptr<Mesh> m_;
	FieldToRCSVector fieldsRCS_;
	std::string basePath_;
	std::vector<std::filesystem::path> dataPaths_;
	double frequency_;

};

}
