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
	/** @brief Performs a RCS calculation of previously exported data through a farfield type probe. The calculated files can be found in data_path in the folders farfield and rcs, with the following naming convention:
	 *  Data_dimensionstring_Th_thetaangle_Phi_phiangle_dgtd.dat
	 * @param[in] data_path Root folder where the farfield probe has exported the simulation data.
	 * @param[in] json_path Root folder with the .json file used to define the simulation.
	 * @param[in] frequencies Standard Library vector of doubles with frequencies defined in the International System (Hz).
	 * @param[in] angles Standard Library vector of SphericalAngles, a struct with theta and phi angles at the desired observation point.
	  */
	RCSManager(const std::string& data_path, const std::string& json_path, std::vector<Frequency>& frequencies, const std::vector<SphericalAngles>& angles);

private:

	std::map<SphericalAngles, Freq2Value> RCSdata_;

};

}
