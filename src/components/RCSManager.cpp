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
		res[f] = physicalConstants::vacuumPermittivity * std::pow(std::abs(freq_val), 2.0);
	}

	//normalize by max val
	auto max = *std::max_element(std::begin(res), std::end(res));
	for (auto f{ 0 }; f < res.size(); f++) {
		res[f] /= max;
	}

	return res;
}


RCSManager::RCSManager(const std::string& data_path, const std::string& json_path, std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{
	auto dim{ Mesh::LoadFromFile(data_path + "/mesh", 1, 0).Dimension() };
	std::string dim_str;
	dim == 2 ? dim_str = "2D_" : dim_str = "3D_";

	std::vector<double> rescaled_frequencies(frequencies.size());
	for (auto f{0}; f < rescaled_frequencies.size(); f++){
		rescaled_frequencies[f] = frequencies[f] / physicalConstants::speedOfLight_SI;
	}

	auto pot_inc{ buildNormalizationTerm(json_path, data_path, rescaled_frequencies) };
	FarField ff(data_path, json_path, rescaled_frequencies, angle_vec);

	if (!std::filesystem::is_directory(data_path + "/farfield") || !std::filesystem::exists(data_path + "/farfield")) {
		std::filesystem::create_directory(data_path + "/farfield");
	}

	if (!std::filesystem::is_directory(data_path + "/rcs") || !std::filesystem::exists(data_path + "/rcs")) {
		std::filesystem::create_directory(data_path + "/rcs");
	}

	double freqdata, const_term, landa, wavenumber;
	for (const auto& angpair : angle_vec) {
		std::ofstream myfile;
		std::string path(data_path + "/farfield/farfieldData_" + dim_str + "Th_" + std::to_string(angpair.theta) + "_" + "Phi_" + std::to_string(angpair.phi) + "_dgtd.dat");
		myfile.open(path);
		if (myfile.is_open()) {
			myfile << "Theta (rad) // " << "Phi (rad) // " << "Frequency (Hz) // " << "r2 * pot // " << "normalization_term\n";
			for (const auto& f : rescaled_frequencies) {
				landa = physicalConstants::speedOfLight / f;
				double normalization;
				dim == 2 ? normalization = landa : normalization = landa * landa;
				myfile << angpair.theta << " " << angpair.phi << " " << f * physicalConstants::speedOfLight_SI << " " << ff.getPotRad(angpair, f) << " " << normalization << "\n";
			}
			myfile.close();
		}
		else {
			throw std::runtime_error("Could not open file to write FarField data.");
		}
	}

	for (int f{ 0 }; f < rescaled_frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			const_term = 4.0 * M_PI / pot_inc[f]; // We defined FarField as r^2 * magnitude, thus we do not include the r^2 found in Equation 8.36 in Taflove's book.
			RCSdata_[angpair][rescaled_frequencies[f]] = const_term * ff.getPotRad(angpair, rescaled_frequencies[f]);
		}
	}

	for (const auto& angpair : angle_vec) {
		std::ofstream myfile;
		std::string path(data_path + "/rcs/rcsData_" + dim_str + "Th_" + std::to_string(angpair.theta) + "_" + "Phi_" + std::to_string(angpair.phi) + "_dgtd.dat");
		myfile.open(path);
		if (myfile.is_open()) {
			myfile <<  "Theta (rad) // " << "Phi (rad) // " << "Frequency (Hz) // " << "rcs // " << "normalization_term\n";
			for (const auto& f : rescaled_frequencies) {
				landa = physicalConstants::speedOfLight / f;
				double normalization;
				dim == 2 ? normalization = landa : normalization = landa * landa; // We will add 2.0 * M_PI to the RCSdata due to the limits regarding a DFT, as we need to increase the range.
				myfile << angpair.theta << " " << angpair.phi << " " << f * physicalConstants::speedOfLight_SI << " " << 2.0 * M_PI * RCSdata_[angpair][f] << " " << normalization << "\n";
			}
			myfile.close();
		}
		else {
			throw std::runtime_error("Could not open file to write RCS data.");

		}
	}

}


}
