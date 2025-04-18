#include "components/RCSManager.h"
#include "math/PhysicalConstants.h"

namespace maxwell {

using namespace mfem;

using json = nlohmann::json;

void FreqFields::append(ComplexVector vector, const std::string& field, const size_t freq)
{
	if (field == "/Ex.gf") {
		Ex[freq] = vector;
	}
	else if (field == "/Ey.gf")
	{
		Ey[freq] = vector;
	}
	else if (field == "/Ez.gf")
	{
		Ez[freq] = vector;
	}
	else if (field == "/Hx.gf")
	{
		Hx[freq] = vector;
	}
	else if (field == "/Hy.gf")
	{
		Hy[freq] = vector;
	}
	else if (field == "/Hz.gf")
	{
		Hz[freq] = vector;
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
	auto landa = physicalConstants::speedOfLight_SI / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * cos(phi) + x[1] * sin(phi));
	return cos(rad_term);
}
double func_exp_imag_part_2D(const Vector& x, const double freq, const Phi phi)
{
	auto landa = physicalConstants::speedOfLight_SI / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * cos(phi) + x[1] * sin(phi));
	return sin(rad_term);
}
double func_exp_real_part_3D(const Vector& x, const double freq, const SphericalAngles angles)
{
	auto landa = physicalConstants::speedOfLight_SI / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto rad_term = wavenumber * (x[0] * sin(angles.second) * cos(angles.first) + x[1] * sin(angles.second) * sin(angles.first) + x[2] * cos(angles.second));
	return cos(rad_term);
}
double func_exp_imag_part_3D(const Vector& x, const double freq, const SphericalAngles angles)
{
	auto landa = physicalConstants::speedOfLight_SI / freq;
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
std::unique_ptr<LinearForm> assembleLinearForm(FunctionCoefficient* fc, FiniteElementSpace& fes, const Direction& dir)
{
	auto res{ std::make_unique<LinearForm>(&fes) };
	auto marker{ getNTFFMarker(fes.GetMesh()->bdr_attributes.Max()) };
	res->AddBdrFaceIntegrator(new mfemExtension::RCSBdrFaceIntegrator(*fc, dir), marker);
	res->Assemble();
	return res;
}
std::unique_ptr<FunctionCoefficient> buildFC_2D(const double freq, const Phi& phi, bool isReal)
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
	return std::make_unique<FunctionCoefficient>(res);
}
std::unique_ptr<FunctionCoefficient> buildFC_3D(const double freq, const SphericalAngles& angles, bool isReal)
{
	std::function<double(const Vector&)> f = 0;
	switch (isReal) {
	case true:
		f = std::bind(&func_exp_real_part_3D, std::placeholders::_1, freq, angles);
		break;
	case false:
		f = std::bind(&func_exp_imag_part_3D, std::placeholders::_1, freq, angles);
		break;
	}
	FunctionCoefficient res(f);
	return std::make_unique<FunctionCoefficient>(res);
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

std::map<SphericalAngles, Freq2Value> initAngles2FreqValues(const std::vector<double>& frequencies, const std::vector<SphericalAngles>& angleVec)
{
	std::map<SphericalAngles, Freq2Value> res;
	for (const auto& angpair : angleVec) {
		Freq2Value f2v;
		for (const auto& f : frequencies) {
			f2v.emplace(f, 0.0);
		}
		res.emplace(angpair, f2v);
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

	if (std::abs(mean - 1e5) > 1e-6 || std::abs(delay - 1e5) > 1e-6) {
		throw std::runtime_error("Verify PlaneWaveData inputs for RCS normalization term.");
	}

	return PlaneWaveData(mean / physicalConstants::speedOfLight_SI, delay / physicalConstants::speedOfLight_SI);
}
std::vector<double> buildTimeVector(const std::string& data_path) 
{
	std::vector<double> res;
	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			res.push_back(getTime(dir_entry.path().generic_string() + "/time.txt") / physicalConstants::speedOfLight_SI);
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
void exportIncidentGaussData(std::vector<double>& time, std::vector<double>& gauss_val, const std::string& json_data)
{
	std::ofstream myfile;
	myfile.open("../personal-sandbox/Python/GaussData_" + json_data + "_dgtd.dat");
	myfile << "Time (s) " << "Gaussian Val (mag)\n";
	for (int t = 0; t < time.size(); ++t) {
		myfile << time[t] << " " << gauss_val[t] << "\n";
	}
	myfile.close();
}
void exportTransformGaussData(std::vector<double>& frequencies, std::map<double, std::complex<double>>& map, const std::string& json_data)
{
	std::ofstream myfile;
	myfile.open("../personal-sandbox/Python/TransformData_" + json_data + "_dgtd.dat");
	myfile << "Frequency (Hz) " << "Real " << "Imag\n";
	for (int f = 0; f < frequencies.size(); ++f) {
		myfile << frequencies[f] << " " << map[frequencies[f]].real() << " " << map[frequencies[f]].imag() << "\n";
	}
	myfile.close();
}

std::vector<double> buildNormalizationTerm(const std::string& json_path, const std::string& path, std::vector<double>& frequencies)
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
	auto max = *std::max_element(std::begin(res), std::end(res));
	for (auto f{ 0 }; f < res.size(); f++) {
		res[f] /= max;
	}

	trimLowMagFreqs(freq2complex, frequencies);
	exportIncidentGaussData(time, gauss_val, json_path);
	exportTransformGaussData(frequencies, freq2complex, json_path);

	return res;
}

Freq2CompVec calculateDFT(const Vector& gf, const std::vector<double>& frequencies, const double time)
{
	Freq2CompVec res(frequencies.size());
	for (int f{ 0 }; f < frequencies.size(); f++) {
		res[f].resize(gf.Size(), std::complex<double>(0.0, 0.0));
		for (int i{ 0 }; i < gf.Size(); i++) {
			auto arg = 2.0 * M_PI * frequencies[f] * time;
			auto w = std::complex<double>(cos(arg), -sin(arg));
			res[f][i] += gf[i] * w;
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
	FreqFields res(frequencies.size());
	std::vector<std::string> fields({ "/Ex.gf", "/Ey.gf", "/Ez.gf", "/Hx.gf", "/Hy.gf", "/Hz.gf" });

	std::vector<double> time{ buildTimeVector(path) };

	for (const auto& field : fields) {
		std::vector<GridFunction> A;
		for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
			if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
				A.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + field));
			}
		}
		if (field == "/Hx.gf" || field == "/Hy.gf" || field == "/Hz.gf") {
			for (int g{ 0 }; g < A.size(); ++g) {
				A[g] /= physicalConstants::freeSpaceImpedance_SI;
			}
		}
		for (int f{ 0 }; f < frequencies.size(); ++f) {
			ComplexVector comp_vec(A[f].Size());
			for (int t{ 0 }; t < time.size(); ++t) {
				auto arg = 2.0 * M_PI * frequencies[f] * time[t];
				auto w = std::complex<double>(cos(arg), -sin(arg));
				for (int v{ 0 }; v < A[f].Size(); ++v) {
					comp_vec[v] += A[t][v] * w;
				}
			}
			res.append(comp_vec, field, f);
		}
	}

	normaliseFreqFields(res, time.size());
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

std::pair<std::complex<double>, std::complex<double>> RCSManager::performRCSCalculations(ComplexVector& FAx, ComplexVector& FAy, ComplexVector& FAz, const double frequency, const SphericalAngles& angles, bool isElectric)
{
	std::unique_ptr<FunctionCoefficient> fc_real, fc_imag;
	
	switch (fes_->GetMesh()->SpaceDimension()) {
	case 2:
		fc_real = buildFC_2D(frequency, angles.first, true) ;
		fc_imag = buildFC_2D(frequency, angles.first, false);
		break;
	case 3:
		fc_real = buildFC_3D(frequency, angles, true) ;
		fc_imag = buildFC_3D(frequency, angles, false);
	}
	
	auto funcCoeff {std::make_pair(fc_real.get(), fc_imag.get())};
	auto lf_x{ assembleComplexLinearForm(funcCoeff, *fes_.get(), X) };
	auto lf_y{ assembleComplexLinearForm(funcCoeff, *fes_.get(), Y) };
	auto lf_z{ assembleComplexLinearForm(funcCoeff, *fes_.get(), Z) };

	auto DCx{ complexInnerProduct(lf_y, FAz) - complexInnerProduct(lf_z, FAy) };
	auto DCy{ complexInnerProduct(lf_z, FAx) - complexInnerProduct(lf_x, FAz) };
	auto DCz{ complexInnerProduct(lf_x, FAy) - complexInnerProduct(lf_y, FAx) };

	if (isElectric) {
		DCx *= -1.0;
		DCy *= -1.0;
		DCz *= -1.0;
	}

	auto phi_value   = -DCx * sin(angles.first)                      + DCy * cos(angles.first);
	auto theta_value =  DCx * cos(angles.second) * cos(angles.first) + DCy * cos(angles.second) * sin(angles.first) - DCz * sin(angles.second);

	return std::pair<std::complex<double>, std::complex<double>>(phi_value, theta_value);
}

RCSManager::RCSManager(const std::string& path, const std::string& json_path, std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{

	Mesh mesh{ Mesh::LoadFromFile(path + "/mesh", 1, 0) };
	getFESFromGF(mesh, path);

	auto RCSdata{ initAngles2FreqValues(frequencies, angle_vec) };

	auto normalization_term{ buildNormalizationTerm(json_path, path, frequencies) };
	
	FreqFields FreqFields{ calculateFreqFields(mesh, frequencies, path)};

	double freqdata, const_term, landa, wavenumber;
	std::pair<std::complex<double>, std::complex<double>> N_pair, L_pair;
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			N_pair = performRCSCalculations(FreqFields.Hx[f], FreqFields.Hy[f], FreqFields.Hz[f], frequencies[f], angpair, false);
			L_pair = performRCSCalculations(FreqFields.Ex[f], FreqFields.Ey[f], FreqFields.Ez[f], frequencies[f], angpair, true);
			landa = physicalConstants::speedOfLight_SI / frequencies[f];
			wavenumber = 2.0 * M_PI / landa;
			const_term = std::pow(wavenumber, 2.0) / (8.0 * M_PI * physicalConstants::freeSpaceImpedance_SI * normalization_term[f]);
			freqdata = const_term * (std::pow(std::abs(L_pair.first + physicalConstants::freeSpaceImpedance_SI * N_pair.second), 2.0) + std::pow(std::abs(L_pair.second - physicalConstants::freeSpaceImpedance_SI * N_pair.first), 2.0));
			RCSdata[angpair][frequencies[f]] = freqdata;
		}
	}

	std::string dim;
	fes_->GetMesh()->SpaceDimension() == 2 ? dim = "2D_" : dim = "3D_";

	for (const auto& angpair : angle_vec) {
		std::ofstream myfile;
		myfile.open("../personal-sandbox/Python/RCSData_" + json_path + dim + std::to_string(angpair.first) + "_" + std::to_string(angpair.second) + "_dgtd.dat");
		myfile << "Angle Rho " << "Angle Phi " << "Frequency (Hz) " << "rcs\n";
		for (const auto& f : frequencies) {
			auto landa = physicalConstants::speedOfLight_SI / f;
			double normalization;
			fes_->GetMesh()->SpaceDimension() == 2 ? normalization = landa : normalization = landa * landa;
			myfile << angpair.first << " " << angpair.second << " " << f << " " << RCSdata[angpair][f] << "\n";
		}
		myfile.close();
	}

}


}
