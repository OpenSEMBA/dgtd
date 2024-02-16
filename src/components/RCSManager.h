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
using SphericalAngles = std::pair<Phi, Rho>;
using ComplexVector = std::vector<std::complex<double>>;
using DFTFreqFieldsComp = std::vector<std::vector<std::complex<double>>>;
using DFTFreqFieldsDouble = std::vector<std::vector<double>>;

struct RCSData {
	double RCSvalue;
	double frequency;
	SphericalAngles angles;
	double time;

	RCSData(double val, double freq, SphericalAngles);
};

class RCSManager {
public:

	RCSManager(const std::string& path, const std::vector<double>& frenquency, const std::vector<SphericalAngles>& angle);

private:

	std::pair<std::complex<double>, std::complex<double>> performRCS2DCalculations(ComplexVector& FAx, ComplexVector& Ay, ComplexVector& Az, const double frequency, const SphericalAngles&);
	DFTFreqFieldsComp assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& field);
	void fillPostDataMaps(const std::vector<double>& frequencies, const std::vector<SphericalAngles>& angleVec);
	
	Mesh m_;
	std::unique_ptr<DG_FECollection> fec_;
	std::unique_ptr<FiniteElementSpace> fes_;
	std::string base_path_;
	std::map<SphericalAngles, Freq2Value> postdata_;

};

}
