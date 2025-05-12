#include "components/RCSManager.h"

namespace maxwell {

using namespace mfem;

using json = nlohmann::json;

static std::vector<double> buildNormalizationTerm(const std::string& json_path, const std::string& path, std::vector<double>& frequencies)
{

	auto case_data = driver::parseJSONfile(json_path);
	double spread = case_data["sources"][0]["magnitude"]["spread"];
	mfem::Vector mean = driver::assemble3DVector(case_data["sources"][0]["magnitude"]["mean"]);
	mfem::Vector propagation = driver::assemble3DVector(case_data["sources"][0]["propagation"]);
	double projMean = mean * propagation / propagation.Norml2();;
	std::vector<double> time{ buildTimeVector(path) };
	std::vector<double> gauss_val{ evaluateGaussianVector(time, spread, projMean) };

	std::map<double, std::complex<double>> freq2complex;
	std::vector<double> res(frequencies.size(), 0.0);
	for (int f{ 0 }; f < frequencies.size(); f++) {
		std::complex<double> freq_val(0.0, 0.0);
		for (int t{ 0 }; t < time.size(); t++) {
			auto arg = 2.0 * M_PI * frequencies[f] * time[t];
			auto w = std::complex<double>(cos(arg), -sin(arg));
			freq_val += gauss_val[t] * w; 
		}
		freq2complex.emplace(std::make_pair(frequencies[f],freq_val));
		res[f] = std::pow(std::abs(freq_val), 2.0) / (2.0 * physicalConstants::freeSpaceImpedance);
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

	std::filesystem::path rcs_path(data_path + "/rcs");
	if (std::filesystem::exists(rcs_path)) {
		std::error_code ec;
		std::filesystem::remove_all(rcs_path, ec);
		if (ec) {
			std::cerr << "Error removing directory: " << ec.message() << '\n';
		}
	}

	FarField ff(data_path, json_path, rescaled_frequencies, angle_vec);
	auto pot_inc{ buildNormalizationTerm(json_path, data_path, rescaled_frequencies) };


	if (!std::filesystem::exists(rcs_path)) {
		std::filesystem::create_directory(rcs_path);
	}

	double landa;

	for (const auto& angpair : angle_vec) {
		std::ofstream myfile;
		std::string path(data_path + "/rcs/rcsData_" + dim_str + "Th_" + std::to_string(angpair.theta) + "_" + "Phi_" + std::to_string(angpair.phi) + "_dgtd.dat");
		myfile.open(path);
		if (myfile.is_open()) {
			myfile <<  "Theta (rad) // " << "Phi (rad) // " << "Frequency (Hz) // " << "rcs // " << "normalization_term\n";
			for (auto f = 0; f < frequencies.size(); f++) {
				landa = physicalConstants::speedOfLight / rescaled_frequencies[f];
				double normalization;
				dim == 2 ? normalization = landa : normalization = landa * landa;
				myfile << angpair.theta << " " << angpair.phi << " " << frequencies[f] << " " << ff.getPotRad(angpair, rescaled_frequencies[f]) / pot_inc[f] << " " << normalization << "\n";
			}
			myfile.close();
		}
		else {
			throw std::runtime_error("Could not open file to write RCS data.");

		}
	}

}


}
