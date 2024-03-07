#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;

using json = nlohmann::json;


//double mu_0 = 1.0;
//double e_0 = 1.0;
double mu_0 = 4.0e-7 * M_PI;
double e_0 = 8.8541878128e-12;
double eta_0{ sqrt(mu_0 / e_0) };
double speed_of_wave{ 1.0 / sqrt(mu_0 * e_0) };

void FreqFields::append(Freq2CompVec f2cv, const std::string& field)
{
	if (field == "/Ex.gf") {
		Ex = f2cv;
	}
	else if (field == "/Ey.gf")
	{
		Ey = f2cv;
	}
	else if (field == "/Ez.gf")
	{
		Ez = f2cv;
	}
	else if (field == "/Hx.gf")
	{
		Hx = f2cv;
	}
	else if (field == "/Hy.gf")
	{
		Hy = f2cv;
	}
	else if (field == "/Hz.gf")
	{
		Hz = f2cv;
	}
};

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
	auto landa = speed_of_wave / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * cos(phi) + x[1] * sin(phi));
	return cos(rad_term);
}
double func_exp_imag_part_2D(const Vector& x, const double freq, const Phi phi)
{
	auto landa = speed_of_wave / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * cos(phi) + x[1] * sin(phi));
	return sin(rad_term);
}
double func_exp_real_part_3D(const Vector& x, const double freq, const SphericalAngles angles)
{
	auto landa = speed_of_wave / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * sin(angles.second) * cos(angles.first) + x[1] * sin(angles.second) * sin(angles.first) + x[2] * cos(angles.second));
	return cos(rad_term);
}
double func_exp_imag_part_3D(const Vector& x, const double freq, const SphericalAngles angles)
{
	auto landa = speed_of_wave / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * sin(angles.second) * cos(angles.first) + x[1] * sin(angles.second) * sin(angles.first) + x[2] * cos(angles.second));
	return sin(rad_term);
}

void RCSManager::getFESFromGF(Mesh& mesh, const std::string& path)
{
	for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			std::ifstream inex(dir_entry.path().generic_string() + "/Ex.gf");
			FiniteElementSpace fes;
			fes.Load(&mesh, inex);
			fes_ = std::make_unique<FiniteElementSpace>(fes);
			break;
		}
	}
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

std::map<SphericalAngles, Freq2Value> fillPostDataMaps(const std::vector<double>& frequencies, const std::vector<SphericalAngles>& angleVec)
{
	std::map<SphericalAngles, Freq2Value> res;
	for (const auto& angpair : angleVec) {
		Freq2Value inner;
		for (const auto& f : frequencies) {
			inner.emplace(f, 0.0);
		}
		res.emplace(angpair, inner);
	}
	return res;
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

	return PlaneWaveData(mean / speed_of_wave, delay / speed_of_wave);
}
std::vector<double> buildTimeVector(const std::string& data_path) 
{
	std::vector<double> res;
	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			res.push_back(getTime(dir_entry.path().generic_string() + "/time.txt") / speed_of_wave);
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

std::vector<double> buildNormalizationTerm(const std::string& json_path, const std::string& path, std::vector<double>& frequencies)
{

	auto planewave_data{ buildPlaneWaveData(parseJSONfile(json_path)) };
	std::vector<double> time{ buildTimeVector(path) };
	std::vector<double> gauss_val{ evaluateGaussianVector(time, planewave_data.delay, planewave_data.mean) };

	std::map<double, std::complex<double>> freq2complex;
	std::vector<double> res(frequencies.size(), 0.0);
	for (int f{ 0 }; f < frequencies.size(); f++) {
		std::complex<double> freq_val(0.0, 0.0);
		for (int t{ 0 }; t < time.size(); t++) {
			auto arg = -2.0 * M_PI * frequencies[f] * time[t];
			freq_val += std::complex<double>(gauss_val[t] * cos(arg), gauss_val[t] * sin(arg));
		}
		freq2complex.emplace(std::make_pair(frequencies[f], freq_val));
		res[f] = e_0 * std::pow(std::abs(freq_val), 2.0);
	}

	trimLowMagFreqs(freq2complex, frequencies);
	exportIncidentGaussData(time, gauss_val);
	exportTransformGaussData(frequencies, freq2complex);

	return res;
}

Freq2CompVec calculateDFT(const Vector& gf, const std::vector<double>& frequencies, const double time)
{
	Freq2CompVec res(frequencies.size());
	for (int f{ 0 }; f < frequencies.size(); f++) {
		res[f].resize(gf.Size());
		for (int i{ 0 }; i < gf.Size(); i++) {
			auto arg = -2.0 * M_PI * frequencies[f] * time;
			res[f][i] += std::complex<double>(gf[i] * cos(arg), gf[i] * sin(arg));
		}
	}
	return res;
}
void normaliseFreqFields(FreqFields& ff, size_t value)
{
	for (int f{ 0 }; f < ff.Ex.size(); ++f) {
		for (int v{ 0 }; v < ff.Ex[f].size(); ++v) {
			ff.Ex[f][v] /= (double)value;
			ff.Ey[f][v] /= (double)value;
			ff.Ez[f][v] /= (double)value;
			ff.Hx[f][v] /= (double)value;
			ff.Hy[f][v] /= (double)value;
			ff.Hz[f][v] /= (double)value;
		}
	}
}
FreqFields calculateFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& path)
{
	FreqFields res;
	size_t time_step_counter{ 0 };
	std::vector<std::string> fields({ "/Ex.gf", "/Ey.gf", "/Ez.gf", "/Hx.gf", "/Hy.gf", "/Hz.gf" });

	std::vector<double> time{ buildTimeVector(path) };

	std::vector<GridFunction> Ex;
	std::vector<GridFunction> Ey;
	std::vector<GridFunction> Ez;
	std::vector<GridFunction> Hx;
	std::vector<GridFunction> Hy;
	std::vector<GridFunction> Hz;

	for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			Ex.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + "/Ex.gf"));
			Ey.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + "/Ey.gf"));
			Ez.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + "/Ez.gf"));
			Hx.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + "/Hx.gf"));
			Hy.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + "/Hy.gf"));
			Hz.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + "/Hz.gf"));
		}
	}

	for (int g{ 0 }; g < Hx.size(); ++g) {
		Hx[g] /= eta_0;
		Hy[g] /= eta_0;
		Hz[g] /= eta_0;
	}

	res.Ex.resize(frequencies.size());
	res.Ey.resize(frequencies.size());
	res.Ez.resize(frequencies.size());
	res.Hx.resize(frequencies.size());
	res.Hy.resize(frequencies.size());
	res.Hz.resize(frequencies.size());

	for (int f{ 0 }; f < frequencies.size(); ++f) {
	ComplexVector comp_ex(Ex[f].Size()), comp_ey(Ey[f].Size()), comp_ez(Ez[f].Size()), comp_hx(Hx[f].Size()), comp_hy(Hy[f].Size()), comp_hz(Hz[f].Size());
		for (int t{ 0 }; t < time.size(); ++t) {
			auto arg = -2.0 * M_PI * frequencies[f] * time[t];
			for (int v{ 0 }; v < Ex[f].Size(); ++v) {
				comp_ex[v] += std::complex<double>(Ex[t][v] * cos(arg), Ex[t][v] * sin(arg));
				comp_ey[v] += std::complex<double>(Ey[t][v] * cos(arg), Ey[t][v] * sin(arg));
				comp_ez[v] += std::complex<double>(Ez[t][v] * cos(arg), Ez[t][v] * sin(arg));
				comp_hx[v] += std::complex<double>(Hx[t][v] * cos(arg), Hx[t][v] * sin(arg));
				comp_hy[v] += std::complex<double>(Hy[t][v] * cos(arg), Hy[t][v] * sin(arg));
				comp_hz[v] += std::complex<double>(Hz[t][v] * cos(arg), Hz[t][v] * sin(arg));
			}
			time_step_counter++;
		}
		res.Ex[f] = comp_ex;
		res.Ey[f] = comp_ey;
		res.Ez[f] = comp_ez;
		res.Hx[f] = comp_hx;
		res.Hy[f] = comp_hy;
		res.Hz[f] = comp_hz;
	}

	normaliseFreqFields(res, time_step_counter / frequencies.size());
	return res;
}

ComplexVector assembleComplexLinearForm(FunctionPair& fp, FiniteElementSpace& fes, const Direction& dir) 
{
	ComplexVector res;
	std::unique_ptr<LinearForm> lf_real = assembleLinearForm(fp.first, fes, dir);
	std::unique_ptr<LinearForm> lf_imag = assembleLinearForm(fp.second, fes, dir);
	res.resize(lf_real->Size());
	for (int i{ 0 }; i < res.size(); i++) {
		res[i] = std::complex<double>(lf_real->Elem(i), lf_imag->Elem(i));
	}
	return res;
}

std::pair<std::complex<double>, std::complex<double>> RCSManager::performRCS2DCalculations(ComplexVector& FAx, ComplexVector& FAy, ComplexVector& FAz, const double frequency, const SphericalAngles& angles, bool isElectric)
{
	auto fc_real{ buildFC_2D(frequency, angles.first, true) };
	auto fc_imag{ buildFC_2D(frequency, angles.first, false) };

	auto lf_x{ assembleComplexLinearForm(std::make_pair(fc_real,fc_imag),*fes_.get(),X) };
	auto lf_y{ assembleComplexLinearForm(std::make_pair(fc_real,fc_imag),*fes_.get(),Y) };
	auto lf_z{ assembleComplexLinearForm(std::make_pair(fc_real,fc_imag),*fes_.get(),Z) };

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

	auto phi_value   = -DCx * sin(angles.first)                      + DCy * cos(angles.first);
	auto theta_value =  DCx * cos(angles.second) * cos(angles.first) + DCy * cos(angles.second) * sin(angles.first) - DCz * sin(angles.second);

	return std::pair<std::complex<double>, std::complex<double>>(phi_value, theta_value);
}

RCSManager::RCSManager(const std::string& path, const std::string& json_path, double dt, int steps, const std::vector<SphericalAngles>& angle_vec)
{

	Mesh mesh{ Mesh::LoadFromFile(path + "/mesh", 1, 0) };
	getFESFromGF(mesh, path);

	const double f_max = 2.0 / dt;
	auto frequencies{ logspace(6.0, 8.0, 10) };
	auto RCSdata{ fillPostDataMaps(frequencies, angle_vec) };

	auto normalization_term{ buildNormalizationTerm(json_path, path, frequencies) };
	
	FreqFields FreqFields{ calculateFreqFields(mesh, frequencies, path)};

	double freqdata, const_term, landa, wavenumber;
	std::pair<std::complex<double>, std::complex<double>> N_pair, L_pair;
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			switch (mesh.SpaceDimension()) {
			case 2:
				N_pair = performRCS2DCalculations(FreqFields.Hx[f], FreqFields.Hy[f], FreqFields.Hz[f], frequencies[f], angpair, false);
				L_pair = performRCS2DCalculations(FreqFields.Ex[f], FreqFields.Ey[f], FreqFields.Ez[f], frequencies[f], angpair, true);
				landa = speed_of_wave / frequencies[f];
				wavenumber = 2.0 * M_PI / landa;
				const_term = std::pow(wavenumber, 2.0) / (8.0 * M_PI * eta_0 * normalization_term[f]);
				freqdata = const_term * (std::pow(std::abs(L_pair.first + eta_0 * N_pair.second), 2.0) + std::pow(std::abs(L_pair.second - eta_0 * N_pair.first), 2.0));
				RCSdata[angpair][frequencies[f]] = freqdata;
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
			auto landa = speed_of_wave / f;
			myfile << angpair.first << " " << angpair.second << " " << f << " " << 10.0*log10(RCSdata[angpair][f]/landa) << "\n";
		}
		myfile.close();
	}

}


}