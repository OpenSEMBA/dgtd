#pragma once

#include <mfem.hpp>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/Probes.h>
#include <mfemExtension/LinearIntegrators.h>


namespace maxwell {

using namespace mfem;
using Rho = double;
using Phi = double;
using Frequency = double;
using RCSValue = double;
using Freq2Value = std::map<Frequency, RCSValue>;
using SphericalAngles = std::pair<Rho, Phi>;
using SphericalAnglesVector = std::vector<SphericalAngles>;
using DFTFreqFieldsComp = std::vector<std::vector<std::complex<double>>>;
using DFTFreqFieldsDouble = std::vector<std::vector<double>>;

struct RCSData {
	double RCSvalue;
	double frequency;
	SphericalAngles angles;
	double time;

	RCSData(double val, double freq, SphericalAngles, double t);
};

class RCSManager {
public:

	RCSManager(const std::string& path, const std::vector<double>& frenquency, const SphericalAnglesVector& angle);

private:

	double performRCS2DCalculations(GridFunction& Ax, GridFunction& Ay, GridFunction& Az, const double frequency, const SphericalAngles&);
	DFTFreqFieldsComp assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& path);
	void fillPostDataMaps(const std::vector<double>& frequencies, const SphericalAnglesVector& angleVec);
	
	Mesh m_;
	std::string base_path_;
	std::map<SphericalAngles, Freq2Value> postdata_;

};

}
