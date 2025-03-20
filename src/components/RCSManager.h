#pragma once

#include "driver/driver.h"
#include "FarField.h"
#include <filesystem>
#include <fstream>

namespace maxwell {

using namespace mfem;

Freq2CompVec calculateDFT(const Vector&, const std::vector<double>& frequencies, const double time);

struct RCSData {
	double RCSvalue;
	double frequency;
	SphericalAngles angles;

	RCSData(double val, double freq, SphericalAngles ang) : 
		RCSvalue(val), 
		frequency(freq),
		angles(ang) 
	{};
};

class RCSManager {
public:

	RCSManager(const std::string& data_path, const std::string& json_path, std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle);

private:

	std::map<SphericalAngles, Freq2Value> RCSdata_;

};

}
