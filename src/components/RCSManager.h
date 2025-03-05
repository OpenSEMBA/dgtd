#pragma once

#include "mfemExtension/LinearIntegrators.h"
#include "driver/driver.h"

namespace maxwell {

using namespace mfem;
using Theta = double;
using Phi = double;
using Frequency = double;
using RCSValue = double;
using Freq2Value = std::map<Frequency, RCSValue>;
using SphericalAngles = std::pair<Phi, Theta>;
using ComplexVector = std::vector<std::complex<double>>;
using Freq2CompVec = std::vector<ComplexVector>;
using DFTFreqFieldsDouble = std::vector<std::vector<double>>;
using FunctionPair = std::pair<FunctionCoefficient*, FunctionCoefficient*>;

Freq2CompVec calculateDFT(const Vector&, const std::vector<double>& frequencies, const double time);

struct FreqFields {

	Freq2CompVec Ex;
	Freq2CompVec Ey;
	Freq2CompVec Ez;
	Freq2CompVec Hx;
	Freq2CompVec Hy;
	Freq2CompVec Hz;

	void append(ComplexVector, const std::string& field, const size_t freq);

	FreqFields(const size_t sizes) {
		Ex.resize(sizes);
		Ey.resize(sizes);
		Ez.resize(sizes);
		Hx.resize(sizes);
		Hy.resize(sizes);
		Hz.resize(sizes);
	}
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

	RCSManager(const std::string& data_path, const std::string& json_path, std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle);

private:

	std::pair<std::complex<double>, std::complex<double>> performRCSCalculations(ComplexVector& FAx, ComplexVector& Ay, ComplexVector& Az, const double frequency, const SphericalAngles&, bool isElectric);
	FreqFields assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& field);
	void getFESFromGF(Mesh& mesh, const std::string& path);

	std::unique_ptr<FiniteElementSpace> fes_;

};

}
