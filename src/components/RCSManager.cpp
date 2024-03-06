#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;

using json = nlohmann::json;

double mu_0 = 1.0;
double e_0 = 1.0;
//double mu_0 = 4.0e-7 * M_PI;
//double e_0 = 8.8541878128e-12;
double eta_0{ sqrt(mu_0 / e_0) };
double v{ 1.0 / sqrt(mu_0 * e_0) };

std::vector<double> logspace(double start, double stop, int num, double base = 10.0) 

{
	std::vector<double> res;
	res.reserve(num);

	double step = (stop - start) / (num - 1);

	for (int i = 0; i < num; ++i) {
		res.push_back(std::pow(base, start + i * step));
	}

	return res;
}
std::complex<double> complexInnerProduct(ComplexVector& first, ComplexVector& second)
{
	if (first.size() != second.size()) {
		throw std::runtime_error("Complex Vectors do not have the same sizes.");
	}
	std::complex<double> res(0.0, 0.0);
	for (int i{ 0 }; i < first.size(); i++) {
		res += first[i] * std::conj(second[i]); //<x,y> = sum(x_i * conj(y_i))
	}
	return res;
}

double func_exp_real_part_2D(const Vector& x, const double freq, const Phi phi)
{
	auto landa = v / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * cos(phi) + x[1] * sin(phi));
	return cos(rad_term);
}
double func_exp_imag_part_2D(const Vector& x, const double freq, const Phi phi)
{
	auto landa = v / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * cos(phi) + x[1] * sin(phi));
	return sin(rad_term);
}
double func_exp_real_part_3D(const Vector& x, const double freq, const SphericalAngles angles)
{
	auto landa = v / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * sin(angles.second) * cos(angles.first) + x[1] * sin(angles.second) * sin(angles.first) + x[2] * cos(angles.second));
	return cos(rad_term);
}
double func_exp_imag_part_3D(const Vector& x, const double freq, const SphericalAngles angles)
{
	auto landa = v / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * sin(angles.second) * cos(angles.first) + x[1] * sin(angles.second) * sin(angles.first) + x[2] * cos(angles.second));
	return sin(rad_term);
}

Array<int> getNTFFMarker(const int att_size)
{
	Array<int> res(att_size);
	res = 0;
	res[static_cast<int>(BdrCond::NearToFarField) - 1] = 1;
	return res;
}
std::unique_ptr<LinearForm> assembleLinearForm(FunctionCoefficient& fc, FiniteElementSpace& fes, const Direction& dir)
{
	auto res{ std::make_unique<LinearForm>(&fes) };
	auto marker{ getNTFFMarker(fes.GetMesh()->bdr_attributes.Max()) };
	res->AddBdrFaceIntegrator(new mfemExtension::RCSBdrFaceIntegrator(fc, dir), marker);
	res->Assemble();
	return res;
}
FunctionCoefficient buildFC_2D(const double freq, const Phi& phi, bool isReal)
{
	std::function<double(const Vector&)> f = 0;
	switch (isReal) {
	case true:
		f = std::bind(&func_exp_real_part_2D, std::placeholders::_1, freq, phi);
		break;
	case false:
		f = std::bind(&func_exp_imag_part_2D, std::placeholders::_1, freq, phi);
		break;
	}
	FunctionCoefficient res(f);
	return res;
}
GridFunction getGridFunction(Mesh& mesh, const std::string& path)
{
	std::ifstream in(path);
	GridFunction res(&mesh, in);
	return res;
}
const double getTime(const std::string& timePath)
{
	std::ifstream timeFile(timePath);
	if (!timeFile) {
		throw std::runtime_error("File could not be opened in getTime.");
	}
	std::string timeString;
	std::getline(timeFile, timeString);
	return std::stod(timeString);
}

void RCSManager::getFESFromGF(Mesh& mesh)
{
	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			std::ifstream inex(dir_entry.path().generic_string() + "/Ex.gf");
			FiniteElementSpace fes;
			fes.Load(&mesh, inex);
			fes_ = std::make_unique<FiniteElementSpace>(fes);
			break;
		}
	}
}
void RCSManager::fillPostDataMaps(const std::vector<double>& frequencies, const std::vector<SphericalAngles>& angleVec)
{
	for (const auto& angpair : angleVec) {
		Freq2Value inner;
		for (const auto& f : frequencies) {
			inner.emplace(f, 0.0);
		}
		postdata_.emplace(angpair, inner);
	}
}


PlaneWaveData buildPlaneWaveData(const json& json)
{
	double mean(-1e5), delay(-1e5);

	for (auto s{ 0 }; s < json["sources"].size(); s++) {
		if (json["sources"][s]["type"] == "totalField") {
			mean = json["sources"][s]["magnitude"]["spread"];
			delay = json["sources"][s]["magnitude"]["delay"];
		}
	}

	if (mean == -1e5 || delay == -1e5) {
		throw std::runtime_error("Verify PlaneWaveData inputs for RCS normalization term.");
	}

	return PlaneWaveData(mean / v, delay / v);
}
std::vector<double> buildTimeVector(const std::string& data_path) 
{
	std::vector<double> res;
	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			res.push_back(getTime(dir_entry.path().generic_string() + "/time.txt") / v);
		}
	}
	return res;
}
std::vector<double> evaluateGaussianVector(std::vector<double>& time, double delay, double mean)
{
	std::vector<double> res(time.size());
	for (int t = 0; t < time.size(); ++t) {
		res[t] = exp(-0.5 * std::pow((time[t] - delay) / mean, 2.0));
	}
	return res;
}
void trimLowMagFreqs(const std::map<double, std::complex<double>>& map, std::vector<double>& frequencies)
{
	const double tol = 1e-2;
	for (int f = 0; f < frequencies.size(); ++f)
	{
		if (std::abs(map.at(frequencies[f])) < tol)
		{
			frequencies.erase(frequencies.begin() + f, frequencies.end());
			break;
		}
	}
}
void exportIncidentGaussData(std::vector<double>& time, std::vector<double>& gauss_val)
{
	std::ofstream myfile;
	myfile.open("../personal-sandbox/Python/GaussData_dgtd.dat");
	myfile << "Time (s) " << "Gaussian Val (mag)\n";
	for (int t = 0; t < time.size(); ++t) {
		myfile << time[t] << " " << gauss_val[t] << "\n";
	}
	myfile.close();
}
void exportTransformGaussData(std::vector<double>& frequencies, std::map<double, std::complex<double>>& map)
{
	std::ofstream myfile;
	myfile.open("../personal-sandbox/Python/TransformData_dgtd.dat");
	myfile << "Frequency (Hz) " << "Real " << "Imag\n";
	for (int f = 0; f < frequencies.size(); ++f) {
		myfile << frequencies[f] << " " << map[frequencies[f]].real() << " " << map[frequencies[f]].imag() << "\n";
	}
	myfile.close();
}

std::vector<double> RCSManager::buildNormalizationTerm(const std::string& json_path, std::vector<double>& frequencies)
{

	auto planewave_data{ buildPlaneWaveData(parseJSONfile(json_path)) };
	std::vector<double> time{ buildTimeVector(data_path_) };
	std::vector<double> gauss_val{ evaluateGaussianVector(time, planewave_data.delay, planewave_data.mean) };

	std::map<double, std::complex<double>> freq2complex;
	std::vector<double> res(frequencies.size(), 0.0);
	for (int f{ 0 }; f < frequencies.size(); f++) {
		std::complex<double> freq_val(0.0, 0.0);
		for (int t{ 0 }; t < time.size(); t++) {
			auto arg = -2.0 * M_PI * frequencies[f] * time[t];
			auto transformed_val = std::complex<double>(gauss_val[t] * cos(arg), gauss_val[t] * sin(arg));
			freq_val += transformed_val;
		}
		freq2complex.emplace(std::make_pair(frequencies[f], freq_val));
		res[f] = e_0 * std::pow(std::abs(freq_val), 2.0);
	}

	trimLowMagFreqs(freq2complex, frequencies);
	exportIncidentGaussData(time, gauss_val);
	exportTransformGaussData(frequencies, freq2complex);

	return res;
}

FreqFields RCSManager::assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& field)
{
	FreqFields res(frequencies.size());
	size_t dofs{ 0 };
	auto time_step_counter{ 0 };

	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			auto gf{ getGridFunction(mesh, dir_entry.path().generic_string() + field) };
			if (field == "/Hx.gf" || field == "/Hy.gf" || field == "/Hz.gf") {
				gf /= eta_0;
			}
			auto time = getTime(dir_entry.path().generic_string() + "/time.txt") / v;
			dofs = gf.Size();
			time_step_counter++;
			for (int f{ 0 }; f < frequencies.size(); f++) {
				if (res[f].size() != dofs) { res[f].resize(dofs); }
				for (int i{ 0 }; i < dofs; i++) {
					auto arg = -2.0 * M_PI * frequencies[f] * time;
					res[f][i] += std::complex<double>(gf[i] * cos(arg), gf[i] * sin(arg));
				}
			}
		}
	}
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (int i{ 0 }; i < dofs; i++) {
			res[f][i] /= time_step_counter;
		}
	}
	return res;
}

std::pair<std::complex<double>, std::complex<double>> RCSManager::performRCS2DCalculations(ComplexVector& FAx, ComplexVector& FAy, ComplexVector& FAz, const double frequency, const SphericalAngles& angles, bool isElectric)
{
	auto fc_real{ buildFC_2D(frequency, angles.first, true) };
	auto fc_imag{ buildFC_2D(frequency, angles.first, false) };

	std::unique_ptr<LinearForm> lf_real_x = assembleLinearForm(fc_real, *fes_.get(), X);
	std::unique_ptr<LinearForm> lf_real_y = assembleLinearForm(fc_real, *fes_.get(), Y);
	std::unique_ptr<LinearForm> lf_real_z = assembleLinearForm(fc_real, *fes_.get(), Z);
	std::unique_ptr<LinearForm> lf_imag_x = assembleLinearForm(fc_imag, *fes_.get(), X);
	std::unique_ptr<LinearForm> lf_imag_y = assembleLinearForm(fc_imag, *fes_.get(), Y);
	std::unique_ptr<LinearForm> lf_imag_z = assembleLinearForm(fc_imag, *fes_.get(), Z);

	ComplexVector lf_x(lf_real_x->Size()), lf_y(lf_real_y->Size()), lf_z(lf_real_z->Size());

	for (int i{ 0 }; i < lf_x.size(); i++) {
		lf_x[i] = std::complex<double>(lf_real_x->Elem(i), lf_imag_x->Elem(i));
		lf_y[i] = std::complex<double>(lf_real_y->Elem(i), lf_imag_y->Elem(i));
		lf_z[i] = std::complex<double>(lf_real_z->Elem(i), lf_imag_z->Elem(i));
	}

	if (isElectric) {
		for (int val{ 0 }; val < FAx.size(); ++val) {
			FAx[val] *= -1.0;
			FAy[val] *= -1.0;
			FAz[val] *= -1.0;
		}
	}

	auto DCx{ complexInnerProduct(lf_y, FAz) - complexInnerProduct(lf_z, FAy) };
	auto DCy{ complexInnerProduct(lf_z, FAx) - complexInnerProduct(lf_x, FAz) };
	auto DCz{ complexInnerProduct(lf_x, FAy) - complexInnerProduct(lf_y, FAx) };

	if (fes_.get()->GetMesh()->SpaceDimension() != 3) {
		DCz = 0.0;
	}

	auto phi_value = -DCx * sin(angles.first)                      + DCy * cos(angles.first);
	auto rho_value =  DCx * cos(angles.second) * cos(angles.first) + DCy * cos(angles.second) * sin(angles.first) - DCz * sin(angles.second);

	return std::pair<std::complex<double>, std::complex<double>>(phi_value, rho_value);
}

RCSManager::RCSManager(const std::string& path, const std::string& json_path, double dt, int steps, const std::vector<SphericalAngles>& angle_vec)
{
	data_path_ = path;
	
	Mesh mesh{ Mesh::LoadFromFile(data_path_ + "/mesh", 1, 0) };
	getFESFromGF(mesh);

	const double f_max = 2.0 / dt;
	auto frequencies{ logspace(-2.0, 2.0, 90) };
	fillPostDataMaps(frequencies, angle_vec);

	auto normalization_term{ buildNormalizationTerm(json_path, frequencies) };
	
	FreqFields FEx{ assembleFreqFields(mesh, frequencies, "/Ex.gf") };
	FreqFields FEy{ assembleFreqFields(mesh, frequencies, "/Ey.gf") };
	FreqFields FEz{ assembleFreqFields(mesh, frequencies, "/Ez.gf") };
 	FreqFields FHx{ assembleFreqFields(mesh, frequencies, "/Hx.gf") };
	FreqFields FHy{ assembleFreqFields(mesh, frequencies, "/Hy.gf") };
	FreqFields FHz{ assembleFreqFields(mesh, frequencies, "/Hz.gf") };

	double freqdata, const_term, landa, wavenumber;
	std::pair<std::complex<double>, std::complex<double>> N_pair, L_pair;
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			switch (mesh.SpaceDimension()) {
			case 2:
				N_pair = performRCS2DCalculations(FHx[f], FHy[f], FHz[f], frequencies[f], angpair, false);
				L_pair = performRCS2DCalculations(FEx[f], FEy[f], FEz[f], frequencies[f], angpair, true);
				landa = v / frequencies[f];
				wavenumber = 2.0 * M_PI / landa;
				const_term = std::pow(wavenumber, 2.0) / (8.0 * M_PI * eta_0 * normalization_term[f]);
				freqdata = const_term * (std::pow(std::abs(L_pair.first + eta_0 * N_pair.second), 2.0) + std::pow(std::abs(L_pair.second - eta_0 * N_pair.first), 2.0));
				postdata_[angpair][frequencies[f]] = freqdata;
				break;
			case 3:
				break;
			default:
				throw std::runtime_error("RCSManager cannot be applied on dimensions other than 2 and 3.");
			}
		}
	}

	for (const auto& angpair : angle_vec) {
		std::ofstream myfile;
		myfile.open("../personal-sandbox/Python/RCSData_" + std::to_string(angpair.first) + "_" + std::to_string(angpair.second) + "_dgtd.dat");
		myfile << "Angle Rho " << "Angle Phi " << "Frequency (Hz) " << "10*log10(RCSData/landa)\n";
		for (const auto& f : frequencies) {
			auto landa = v / f;
			myfile << angpair.first << " " << angpair.second << " " << f << " " << 10.0*log10(postdata_[angpair][f]/landa) << "\n";
		}
		myfile.close();
	}

}


}