#pragma once

#include <mfem.hpp>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/Probes.h>
#include <mfemExtension/LinearIntegrators.h>
#include <adapter/MaxwellAdapter.hpp>

#include <nlohmann/json.hpp>

namespace maxwell {

using namespace mfem;
using Rho = double;
using Phi = double;
using Frequency = double;
using RCSValue = double;
using Freq2Value = std::map<Frequency, RCSValue>;
using SphericalAngles = std::pair<Phi, Rho>;
using ComplexVector = std::vector<std::complex<double>>;
using Freq2CompVec = std::vector<ComplexVector>;
using DFTFreqFieldsDouble = std::vector<std::vector<double>>;
using FunctionPair = std::pair<FunctionCoefficient, FunctionCoefficient>;

struct FreqFields {

	Freq2CompVec Ex;
	Freq2CompVec Ey;
	Freq2CompVec Ez;
	Freq2CompVec Hx;
	Freq2CompVec Hy;
	Freq2CompVec Hz;

	void append(Freq2CompVec, const std::string& field);
};

struct PlaneWaveData {
	double mean;
	double delay;

	PlaneWaveData(double m, double dl) :
		mean(m),
		delay(dl) {};
};

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

	RCSManager(const std::string& data_path, const std::string& json_path, double f_max, int steps, const std::vector<SphericalAngles>& angle);

private:

	std::pair<std::complex<double>, std::complex<double>> performRCS2DCalculations(ComplexVector& FAx, ComplexVector& Ay, ComplexVector& Az, const double frequency, const SphericalAngles&, bool isElectric);
	FreqFields assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& field);
	void getFESFromGF(Mesh& mesh, const std::string& path);

	std::unique_ptr<FiniteElementSpace> fes_;

};

}
