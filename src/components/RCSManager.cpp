#include "components/RCSManager.h"
#include <filesystem>
#include <fstream>
#include <regex>

namespace maxwell {

using namespace mfem;
using json = nlohmann::json;

struct RowData {
    double theta, phi, frequency;
    double pot_sum, norm_sum;
};

using FileRows = std::vector<RowData>;

static std::vector<double> buildIncomingPowerTerm(const std::string& json_path, const std::string& path, std::vector<double>& frequencies)
{

	auto case_data = driver::parseJSONfile(json_path);
	double spread = case_data["sources"][0]["magnitude"]["spread"];
	mfem::Vector mean = driver::assemble3DVector(case_data["sources"][0]["magnitude"]["mean"]);
	mfem::Vector propagation = driver::assemble3DVector(case_data["sources"][0]["propagation"]);
	double projMean = mean * propagation / propagation.Norml2();
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
		freq_val /= (double)time.size();
		freq2complex.emplace(std::make_pair(frequencies[f],freq_val));
		res[f] = std::norm(freq_val) / (2.0 * physicalConstants::freeSpaceImpedance);
	}

	return res;
}

void checkAnglesFor2DProblems(const int dim, const std::vector<SphericalAngles>& angle_vec)
{
	if (dim == 2) {
		for (const auto& angles : angle_vec) {
			if (angles.theta != M_PI_2) {
				throw std::runtime_error("For 2D RCS Post-Processing, angle theta must be pi/2 in every angle pair.");
			}
		}
	}
}

void cleanOutputFolder(const std::string& path)
{
	if (std::filesystem::is_directory(path + "/farfield") || std::filesystem::exists(path + "/farfield")) {
		std::filesystem::remove_all(path + "/farfield");
	}

	if (std::filesystem::is_directory(path + "/rcs") || std::filesystem::exists(path + "/rcs")) {
		std::filesystem::remove_all(path + "/rcs");
	}
}

void createOutputFolder(const std::string& path)
{
	if (!std::filesystem::is_directory(path + "/farfield") || !std::filesystem::exists(path + "/farfield")) {
		std::filesystem::create_directory(path + "/farfield");
	}

	if (!std::filesystem::is_directory(path + "/rcs") || !std::filesystem::exists(path + "/rcs")) {
		std::filesystem::create_directory(path + "/rcs");
	}
}

std::vector<std::string> findRankFolders(const std::string& data_path) {
    std::vector<std::string> rank_paths;

    const std::filesystem::path base_path(data_path);
    if (!std::filesystem::exists(base_path) || !std::filesystem::is_directory(base_path)) {
        std::cerr << "Invalid base path: " << data_path << '\n';
        return rank_paths;
    }

    std::regex rank_dir_pattern(R"(rank\d+)");
    for (const auto& entry : std::filesystem::directory_iterator(base_path)) {
        if (entry.is_directory()) {
            std::string dirname = entry.path().filename().string();
            if (std::regex_match(dirname, rank_dir_pattern)) {
                rank_paths.push_back(entry.path().string());
            }
        }
    }

    return rank_paths;
}

FileRows readDataFile(const std::filesystem::path& filepath) {
    FileRows rows;
    std::ifstream in(filepath);
    if (!in) throw std::runtime_error("Failed to open " + filepath.string());

	std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("Empty file: " + filepath.string()); //Checks for emptiness of file and also helps skip first line which is the header.
    }

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        RowData row;
        if (!(iss >> row.theta >> row.phi >> row.frequency >> row.pot_sum >> row.norm_sum)) {
            throw std::runtime_error("Malformed line in " + filepath.string());
        }
        rows.push_back(row);
    }
    return rows;
}

void addRows(FileRows& base, const FileRows& additional) {
    if (base.size() != additional.size())
        throw std::runtime_error("Mismatch in row count across files.");

    for (size_t i = 0; i < base.size(); ++i) {
        base[i].pot_sum += additional[i].pot_sum;
        base[i].norm_sum += additional[i].norm_sum;
    }
}

void writeMergedFile(const FileRows& rows, const std::filesystem::path& output_path) {
    std::ofstream out(output_path);
	out << "Theta (rad) // " << "Phi (rad) // " << "Frequency (Hz) // " << "pot // " << "normalization_term\n";
    for (const auto& row : rows) {
        out << row.theta << " " << row.phi << " " << row.frequency << " "
            << row.pot_sum << " " << row.norm_sum << "\n";
    }
}

void writeParallelDataIntoSingleFiles(const std::string& case_path) {
    const std::filesystem::path base(case_path);
    std::regex rank_pattern(R"(rank\d+)");
    std::vector<std::filesystem::path> rank_dirs;

    // Discover all rank directories
    for (const auto& entry : std::filesystem::directory_iterator(base)) {
        if (entry.is_directory() &&
            std::regex_match(entry.path().filename().string(), rank_pattern)) {
            rank_dirs.push_back(entry.path());
        }
    }

    if (rank_dirs.empty()) {
        throw std::runtime_error("No rank folders found.");
    }

    for (const std::string& type : { "farfield", "rcs" }) {
        std::filesystem::path output_dir = base / type;
        std::filesystem::create_directories(output_dir);

        std::filesystem::path reference_dir = rank_dirs[0] / type;

        for (const auto& file_entry : std::filesystem::directory_iterator(reference_dir)) {
            if (!file_entry.is_regular_file()) continue;

            const std::string filename = file_entry.path().filename().string();

            FileRows merged_rows = readDataFile(file_entry.path());

            for (size_t i = 1; i < rank_dirs.size(); ++i) {
                std::filesystem::path other_file = rank_dirs[i] / type / filename;
                FileRows other_rows = readDataFile(other_file);
                addRows(merged_rows, other_rows);
            }

            writeMergedFile(merged_rows, output_dir / filename);
        }
    }
}

RCSManager::RCSManager(const std::string& data_path, const std::string& json_path, std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{
	cleanOutputFolder(data_path);

	auto rank_paths = findRankFolders(data_path);
	
	for (auto r = 0; r < rank_paths.size(); r++){
		auto dim{ Mesh::LoadFromFile(rank_paths[r] + "/mesh", 1, 0).Dimension() };
		cleanOutputFolder(rank_paths[r]);
		createOutputFolder(rank_paths[r]);
		checkAnglesFor2DProblems(dim, angle_vec);

		std::vector<double> rescaled_frequencies(frequencies.size());
		for (auto f{0}; f < rescaled_frequencies.size(); f++){
			rescaled_frequencies[f] = frequencies[f] / physicalConstants::speedOfLight_SI;
		}

		FarField ff(rank_paths[r], json_path, rescaled_frequencies, angle_vec);
		auto pot_inc{ buildIncomingPowerTerm(json_path, rank_paths[r], rescaled_frequencies) };

		double landa;
		for (const auto& angpair : angle_vec) {
			std::ofstream myfile;
			std::string path(rank_paths[r] + "/farfield/farfieldData_Th_" + std::to_string(angpair.theta) + "_" + "Phi_" + std::to_string(angpair.phi) + "_dgtd.dat");
			myfile.open(path);
			if (myfile.is_open()) {
				myfile << "Theta (rad) // " << "Phi (rad) // " << "Frequency (Hz) // " << "pot // " << "normalization_term\n";
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
				auto const_term = 4.0 * M_PI / pot_inc[f];
				// dim == 2 ? const_term /= 2.0 : const_term /= M_PI; // If problem is 2D, we divide by the visible length of the circle area of rad 1 (2.0). 
																	// If 3D, by the equivalent surface area of the visible face of the sphere of rad 1 (PI).
				RCSdata_[angpair][rescaled_frequencies[f]] = const_term * ff.getPotRad(angpair, rescaled_frequencies[f]);
			}
		}

		for (const auto& angpair : angle_vec) {
			std::ofstream myfile;
			std::string path(rank_paths[r] + "/rcs/rcsData_Th_" + std::to_string(angpair.theta) + "_" + "Phi_" + std::to_string(angpair.phi) + "_dgtd.dat");
			myfile.open(path);
			if (myfile.is_open()) {
				myfile <<  "Theta (rad) // " << "Phi (rad) // " << "Frequency (Hz) // " << "rcs // " << "normalization_term\n";
				for (const auto& f : rescaled_frequencies) {
					landa = physicalConstants::speedOfLight / f;
					double normalization;
					dim == 2 ? normalization = landa : normalization = landa * landa;
					myfile << angpair.theta << " " << angpair.phi << " " << f * physicalConstants::speedOfLight_SI << " " << RCSdata_[angpair][f] << " " << normalization << "\n";
				}
				myfile.close();
			}
			else {
				throw std::runtime_error("Could not open file to write RCS data.");

			}
		}
	}

	createOutputFolder(data_path);
	writeParallelDataIntoSingleFiles(data_path);
	for (auto r = 0; r < rank_paths.size(); r++){
		cleanOutputFolder(rank_paths[r]);
	}

}


}
