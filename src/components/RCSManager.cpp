#include "components/RCSManager.h"

namespace maxwell {

using namespace mfem;

using json = nlohmann::json;

static std::vector<double> buildNormalizationTerm(const std::string& json_path, const std::string& path, std::vector<double>& frequencies)
{

	auto planewave_data{ buildPlaneWaveData(driver::parseJSONfile(json_path)) };
	std::vector<double> time{ buildTimeVector(path) };
	std::vector<double> gauss_val{ evaluateGaussianVector(time, planewave_data.delay, planewave_data.mean) };

	std::map<double, std::complex<double>> freq2complex;
	std::vector<double> res(frequencies.size(), 0.0);
	for (int f{ 0 }; f < frequencies.size(); f++) {
		std::complex<double> freq_val(0.0, 0.0);
		for (int t{ 0 }; t < time.size(); t++) {
			auto arg = 2.0 * M_PI * frequencies[f] * time[t];
			auto w = std::complex<double>(cos(arg), -sin(arg));
			freq_val += gauss_val[t] * w;
		}
		freq2complex.emplace(std::make_pair(frequencies[f], freq_val));
		res[f] = physicalConstants::vacuumPermittivity_SI * std::pow(std::abs(freq_val), 2.0);
	}

	//normalize by max val
	auto max = *std::max_element(std::begin(res), std::end(res));
	for (auto f{ 0 }; f < res.size(); f++) {
		res[f] /= max;
	}

	//trimLowMagFreqs(freq2complex, frequencies);
	//exportIncidentGaussData(time, gauss_val, json_path);
	//exportTransformGaussData(frequencies, freq2complex, json_path);

	return res;
}


RCSManager::RCSManager(const std::string& data_path, const std::string& json_path, std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{

	auto pot_inc{ buildNormalizationTerm(json_path, data_path, frequencies) };

	FarField ff(data_path, json_path, frequencies, angle_vec);

	double freqdata, const_term, landa, wavenumber;
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			const_term = 4.0 * M_PI / pot_inc[f];
			RCSdata_[angpair][frequencies[f]] = const_term * ff.getPotRad(angpair, frequencies[f]);
		}
	}

	//std::string dim;
	//fes_->GetMesh()->SpaceDimension() == 2 ? dim = "2D_" : dim = "3D_";

	//for (const auto& angpair : angle_vec) {
	//	std::ofstream myfile;
	//	myfile.open("../personal-sandbox/Python/RCSData_" + json_path + dim + std::to_string(angpair.first) + "_" + std::to_string(angpair.second) + "_dgtd.dat");
	//	myfile << "Angle Rho " << "Angle Phi " << "Frequency (Hz) " << "rcs\n";
	//	for (const auto& f : frequencies) {
	//		auto landa = physicalConstants::speedOfLight_SI / f;
	//		double normalization;
	//		fes_->GetMesh()->SpaceDimension() == 2 ? normalization = landa : normalization = landa * landa;
	//		myfile << angpair.first << " " << angpair.second << " " << f << " " << RCSdata[angpair][f] << "\n";
	//	}
	//	myfile.close();
	//}

}


}
