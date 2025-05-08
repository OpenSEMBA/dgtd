#pragma once

#include "mfem.hpp"

#include "components/Types.h"

#include "driver/driver.h"

#include "math/PhysicalConstants.h"

#include "mfemExtension/LinearIntegrators.h"

#include "components/SubMesher.h"
#include "solver/ProbesManager.h"

namespace maxwell {

using namespace mfem;

using Theta = double;
using Phi = double;
using Frequency = double;
using ScPot = double;
using Freq2Value = std::map<Frequency, ScPot>;
using ComplexVector = std::vector<std::complex<double>>;
using Freq2CompVec = std::vector<ComplexVector>;
using DFTFreqFieldsDouble = std::vector<std::vector<double>>;
using FunctionPair = std::pair<VectorFunctionCoefficient*, VectorFunctionCoefficient*>;

const auto obs_radius{ 1e6 }; //Observation point radius distance. We consider a fairly distant enough radius.

struct SphericalAngles {
	double theta;
	double phi;

	bool operator < (const SphericalAngles& angles) const {
		return std::tie(theta, phi) < std::tie(angles.theta, angles.phi);
	}
};

struct PlaneWaveData {
	double spread;
	double mean;

	PlaneWaveData(double s, double m) :
		spread(s),
		mean(m) {};
};

struct NedelecFields {
	std::vector<GridFunction> electric;
	std::vector<GridFunction> magnetic;

	NedelecFields(int size) {
		electric.reserve(size);
		magnetic.reserve(size);
	}
};

struct FreqNedelecComponents {
	GridFunction electric;
	GridFunction magnetic;
};

using RealImagFreqFields = std::pair<FreqNedelecComponents, FreqNedelecComponents>;
using RealImagFreqCurrents = std::pair<FreqNedelecComponents, FreqNedelecComponents>;
using Freq2NedFields = std::map<Frequency, RealImagFreqFields>;

struct ThetaPhiTerms{
	std::complex<double> theta;
	std::complex<double> phi;
};

struct NLTerms {
	ThetaPhiTerms N;
	ThetaPhiTerms L;
};

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

	void normaliseFields(const double val);
};

//double func_exp_real_part_2D(const Position&, const Frequency, const SphericalAngles&);
//double func_exp_imag_part_2D(const Position&, const Frequency, const SphericalAngles&);
void func_exp_real_part_3D(const Position&, Vector& out, const Frequency, const SphericalAngles&);
void func_exp_imag_part_3D(const Position&, Vector& out, const Frequency, const SphericalAngles&);

//std::unique_ptr<FunctionCoefficient> buildFC_2D(const Frequency, const SphericalAngles&, bool isReal);
std::unique_ptr<VectorFunctionCoefficient> buildVFC_3D(const Frequency, const SphericalAngles&, bool isReal);

double calcPsiAngle3D(const Vector& p, const SphericalAngles& angles);

std::complex<double> complexInnerProduct(ComplexVector& first, ComplexVector& second);

std::unique_ptr<FiniteElementSpace> buildFESFromGF(Mesh&, const std::string& data_path);

std::map<SphericalAngles, Freq2Value> initAngles2FreqValues(const std::vector<Frequency>&, const std::vector<SphericalAngles>&);

PlaneWaveData buildPlaneWaveData(const json&);
std::vector<double> buildTimeVector(const std::string& data_path);

GridFunction getGridFunction(Mesh&, const std::string& data_path);
const Time getTime(const std::string& timePath);
std::vector<double> evaluateGaussianVector(std::vector<Time>& time, double delay, double mean);
void trimLowMagFreqs(const std::map<double, std::complex<double>>& map, std::vector<Frequency>&);

Freq2CompVec calculateDFT(const Vector& gf, const std::vector<Frequency>&, const Time);

ComplexVector assembleComplexLinearForm(FunctionPair& fp, FiniteElementSpace&);

Array<int> getNearToFarFieldMarker(const int att_size);

std::unique_ptr<LinearForm> assembleLinearForm(VectorFunctionCoefficient& fc, FiniteElementSpace& fes);

class FarField {
public:
	FarField(const std::string& data_path, const std::string& json_path, std::vector<Frequency>&, const std::vector<SphericalAngles>& angle_vec);
	const double getPotRad(const SphericalAngles& angpair, const Frequency& freq) const { return pot_rad_.at(angpair).at(freq); }

private:

	NedelecFields buildNedelecFields(Mesh& mesh, FiniteElementSpace& dgfes, const std::string& path);
	RealImagFreqCurrents integrateCurrents(const RealImagFreqFields&, const Frequency, const SphericalAngles& angles);
	RealImagFreqFields buildFrequencyFields(const NedelecFields& time_fields, const std::vector<Time>& time, const Frequency frequency);

	std::unique_ptr<FiniteElementSpace> dgfes_, ndfes_;
	std::map<SphericalAngles, Freq2Value> pot_rad_;
};

}