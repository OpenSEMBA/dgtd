#include "components/FarField.h"

namespace maxwell {

using namespace mfem;

const Vector buildObsPointVec(const SphericalAngles& angles) 
{
	const auto r{ obs_radius }; //Observation point radius distance. We consider a fairly distant enough radius.
	const auto& phi{ angles.phi };
	const auto& theta{ angles.theta };
	const auto x{ r * std::sin(theta) * std::cos(phi) };
	const auto y{ r * std::sin(theta) * std::sin(phi) };
	const auto z{ r * std::cos(theta) };
	const Vector res({ x, y, z });
	return res;
}

double calcPsiAngle2D(const Vector& p, const SphericalAngles& angles)
{
	auto vec{ buildObsPointVec(angles) };
	Vector surfVec(vec.Size());
	surfVec[0] = p[0];
	surfVec[1] = p[1];
	surfVec[2] = 0.0;
	return std::acos((vec[0] * surfVec[0] + vec[1] * surfVec[1] + vec[2] * surfVec[2])  / (vec.Norml2() * surfVec.Norml2()));
}

double calcPsiAngle3D(const Vector& p, const SphericalAngles& angles)
{
	auto vec{ buildObsPointVec(angles) };
	return std::acos((vec[0] * p[0] + vec[1] * p[1] + vec[2] * p[2]) / (vec.Norml2() * p.Norml2()));
}

/* These functions represent the exponential part of e^(+j k r' cos(psi)) found in Equations 8.31 and 8.32 of Taflove's Computational Electrodynamics: The Finite-Difference Time-Domain Method. 
   Due to MFEM not supporting complex numbers in its Coefficient class, we separate real and imaginary parts for either 2D or 3D, returning cosine or sine as needed. */
double func_exp_real_part_2D(const Vector& p, const Frequency freq, const SphericalAngles& angles)
{
	auto landa = physicalConstants::speedOfLight / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto psi = calcPsiAngle2D(p, angles);
	auto rad_term = wavenumber * p.Norml2() * std::cos(psi);
	return cos(rad_term);
}

double func_exp_imag_part_2D(const Vector& p, const Frequency freq, const SphericalAngles& angles)
{
	auto landa = physicalConstants::speedOfLight / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto psi = calcPsiAngle2D(p, angles);
	auto rad_term = wavenumber * p.Norml2() * std::cos(psi);
	return sin(rad_term);
}

double func_exp_real_part_3D(const Vector& p, const Frequency freq, const SphericalAngles& angles)
{
	auto landa = physicalConstants::speedOfLight / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto psi = calcPsiAngle3D(p, angles);
	auto rad_term = wavenumber * p.Norml2() * std::cos(psi);
	return cos(rad_term);
}

double func_exp_imag_part_3D(const Vector& p, const Frequency freq, const SphericalAngles& angles)
{
	auto landa = physicalConstants::speedOfLight / freq;
	auto wavenumber = 2.0 * M_PI / landa;
	auto psi = calcPsiAngle3D(p, angles);
	auto rad_term = wavenumber * p.Norml2() * std::cos(psi);
	return sin(rad_term);
}

std::unique_ptr<FunctionCoefficient> buildFC_2D(const Frequency freq, const SphericalAngles& angles, bool isReal)
{
	std::function<double(const Vector&)> f = 0;
	switch (isReal) {
	case true:
		f = std::bind(&func_exp_real_part_2D, std::placeholders::_1, freq, angles);
		break;
	case false:
		f = std::bind(&func_exp_imag_part_2D, std::placeholders::_1, freq, angles);
		break;
	}
	FunctionCoefficient res(f);
	return std::make_unique<FunctionCoefficient>(res);
}

std::unique_ptr<FunctionCoefficient> buildFC_3D(const Frequency freq, const SphericalAngles& angles, bool isReal)
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

std::complex<double> complexInnerProduct(ComplexVector& first, ComplexVector& second)
{
	if (first.size() != second.size()) {
		throw std::runtime_error("Complex Vectors do not have the same sizes.");
	}
	std::complex<double> res(0.0, 0.0);
	for (int i{ 0 }; i < first.size(); i++) {
		res += first[i] * second[i];
	}
	return res;
}

Array<int> getNearToFarFieldMarker(const int att_size)
{
	Array<int> res(att_size);
	res = 0;
	res[static_cast<int>(BdrCond::NearToFarField) - 1] = 1;
	return res;
}

std::unique_ptr<LinearForm> assembleLinearForm(FunctionCoefficient& fc, FiniteElementSpace& fes, const Direction& dir)
{
	auto res{ std::make_unique<LinearForm>(&fes) };
	auto marker{ getNearToFarFieldMarker(fes.GetMesh()->bdr_attributes.Max()) };
	Direction final_dir = dir;
	switch (fes.GetMesh()->Dimension()) {
	case 2:
		if (dir == Z) {
			final_dir = X;
			res->AddBdrFaceIntegrator(new mfemExtension::FarFieldBdrFaceIntegrator(fc, final_dir), marker);
			res->Assemble();
			res.get()->Set(0.0, *res.get());
			break;
		}
		else {
			res->AddBdrFaceIntegrator(new mfemExtension::FarFieldBdrFaceIntegrator(fc, final_dir), marker);
			res->Assemble();
			break;
		}
		break;
	case 3:
		res->AddBdrFaceIntegrator(new mfemExtension::FarFieldBdrFaceIntegrator(fc, final_dir), marker);
		res->Assemble();
		break;
	default:
		throw std::runtime_error("RCS Post processing only supported for dimensions 2 and 3.");
	}
	return res;
}

std::unique_ptr<FiniteElementSpace> buildFESFromGF(Mesh& mesh, const std::string& path)
{
	for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			std::ifstream inex(dir_entry.path().generic_string() + "/Ex.gf");
			FiniteElementSpace fes;
			fes.Load(&mesh, inex);
			return std::make_unique<FiniteElementSpace>(fes);
		}
	}
	throw std::runtime_error("FES from GridFunciton not returned. Verify input path.");
}

std::map<SphericalAngles, Freq2Value> initAngles2FreqValues(const std::vector<Frequency>& frequencies, const std::vector<SphericalAngles>& angleVec)
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

GridFunction getGridFunction(Mesh& mesh, const std::string& path)
{
	std::ifstream in(path);
	GridFunction res(&mesh, in);
	return res;
}

const Time getTime(const std::string& timePath)
{
	std::ifstream timeFile(timePath);
	if (!timeFile) {
		throw std::runtime_error("File could not be opened in getTime.");
	}
	std::string timeString;
	std::getline(timeFile, timeString);
	return std::stod(timeString);
}

PlaneWaveData buildPlaneWaveData(const json& json)
{
	double spread(-1e5);
	mfem::Vector mean;
	double projMean(-1e5);

	for (auto s{ 0 }; s < json["sources"].size(); s++) {
		if (json["sources"][s]["type"] == "planewave") {
			spread = json["sources"][s]["magnitude"]["spread"];
			mean = driver::assemble3DVector(json["sources"][s]["magnitude"]["mean"]);
		    mfem::Vector propagation = driver::assemble3DVector(json["sources"][s]["propagation"]);
			projMean = mean * propagation / propagation.Norml2();
		}
	}

	if (std::abs(spread - 1e5) < 1e-6 || std::abs(projMean - 1e5) < 1e-6) {
		throw std::runtime_error("Verify PlaneWaveData inputs for RCS normalization term.");
	}

	return PlaneWaveData(spread, projMean);
}

std::vector<double> buildTimeVector(const std::string& data_path)
{
	std::vector<double> res;
	for (auto const& dir_entry : std::filesystem::directory_iterator(data_path)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh" && 
			dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 3) != "rcs"  &&
			dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 8) != "farfield") 
		{
			res.push_back(getTime(dir_entry.path().generic_string() + "/time.txt") / physicalConstants::speedOfLight);
		}
	}
	return res;
}

std::vector<double> evaluateGaussianVector(std::vector<Time>& time, double spread, double mean)
{
	std::vector<double> res(time.size());
	for (int t = 0; t < time.size(); ++t) {
		res[t] = std::exp(-std::pow(time[t] - std::abs(mean), 2.0) / (2.0 * std::pow(spread, 2.0)));
	}
	return res;
}

void trimLowMagFreqs(const std::map<Frequency, std::complex<double>>& map, std::vector<Frequency>& frequencies)
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

Freq2CompVec calculateDFT(const Vector& gf, const std::vector<Frequency>& frequencies, const Time time)
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

FreqFields calculateFreqFields(Mesh& mesh, const std::vector<Frequency>& frequencies, const std::string& path)
{
	FreqFields res(frequencies.size());
	std::vector<std::string> fields({ "/Ex.gf", "/Ey.gf", "/Ez.gf", "/Hx.gf", "/Hy.gf", "/Hz.gf" });

	std::vector<double> time{ buildTimeVector(path) };

	for (const auto& field : fields) {
		std::vector<GridFunction> A;
		for (auto const& dir_entry : std::filesystem::directory_iterator(path)) {
			if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh" &&
				dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "/rcs" &&
				dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 9) != "/farfield") {
				A.push_back(getGridFunction(mesh, dir_entry.path().generic_string() + field));
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

	res.normaliseFields(time.size());

	return res;
}

ComplexVector assembleComplexLinearForm(FunctionPair& fp, FiniteElementSpace& fes, const Direction& dir)
{
	ComplexVector res;
	std::unique_ptr<LinearForm> lf_real = assembleLinearForm(*fp.first, fes, dir);
	std::unique_ptr<LinearForm> lf_imag = assembleLinearForm(*fp.second, fes, dir);
	res.resize(lf_real->Size());
	for (int i{ 0 }; i < res.size(); i++) {
		res[i] = std::complex<double>(lf_real->Elem(i), lf_imag->Elem(i));
	}
	return res;
}


std::pair<std::complex<double>, std::complex<double>> FarField::calcNLpair(ComplexVector& FAx, ComplexVector& FAy, ComplexVector& FAz, const Frequency frequency, const SphericalAngles& angles, bool isElectric)
{
	std::unique_ptr<FunctionCoefficient> fc_real, fc_imag;

	switch (fes_->GetMesh()->SpaceDimension()) {
	case 2:
		fc_real = buildFC_2D(frequency, angles, true);
		fc_imag = buildFC_2D(frequency, angles, false);
		break;
	case 3:
		fc_real = buildFC_3D(frequency, angles, true);
		fc_imag = buildFC_3D(frequency, angles, false);
	}

	auto funcCoeff{ std::make_pair(fc_real.get(), fc_imag.get())};
	auto lf_x{ assembleComplexLinearForm(funcCoeff, *fes_.get(), X) };
	auto lf_y{ assembleComplexLinearForm(funcCoeff, *fes_.get(), Y) };
	auto lf_z{ assembleComplexLinearForm(funcCoeff, *fes_.get(), Z) };

	auto DCx{ complexInnerProduct(lf_y, FAz) - complexInnerProduct(lf_z, FAy) };
	auto DCy{ complexInnerProduct(lf_z, FAx) - complexInnerProduct(lf_x, FAz) };
	auto DCz{ complexInnerProduct(lf_x, FAy) - complexInnerProduct(lf_y, FAx) };

	// Currents are defined as J = n x H and M = -n x E. We change the sign of M due to the normal here.
	if (isElectric) {
		DCx *= -1.0;
		DCy *= -1.0;
		DCz *= -1.0;
	}

	auto phi_value = -DCx * sin(angles.phi) + DCy * cos(angles.phi);
	auto theta_value = DCx * cos(angles.theta) * cos(angles.phi) + DCy * cos(angles.theta) * sin(angles.phi) - DCz * sin(angles.theta);

	return std::pair<std::complex<double>, std::complex<double>>(theta_value, phi_value);
}

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

void FreqFields::normaliseFields(const double value)
{
	for (int f{ 0 }; f < this->Ex.size(); ++f) {
		for (int v{ 0 }; v < this->Ex[f].size(); ++v) {
			this->Ex[f][v] /= (double)value;
			this->Ey[f][v] /= (double)value;
			this->Ez[f][v] /= (double)value;
			this->Hx[f][v] /= (double)value;
			this->Hy[f][v] /= (double)value;
			this->Hz[f][v] /= (double)value;
		}
	}
}

FarField::FarField(const std::string& data_path, const std::string& json_path, std::vector<Frequency>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{
	auto mesh = Mesh::LoadFromFile(data_path + "/mesh");
	fes_ = buildFESFromGF(mesh, data_path);

	pot_rad_ = initAngles2FreqValues(frequencies, angle_vec);

	FreqFields freqfields{ calculateFreqFields(mesh, frequencies, data_path) };

	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			auto N_pair = calcNLpair(freqfields.Hx[f], freqfields.Hy[f], freqfields.Hz[f], frequencies[f], angpair, false);
			auto L_pair = calcNLpair(freqfields.Ex[f], freqfields.Ey[f], freqfields.Ez[f], frequencies[f], angpair, true);
			auto landa = physicalConstants::speedOfLight / frequencies[f];
			auto wavenumber = 2.0 * M_PI / landa;
			auto const_term = std::pow(wavenumber, 2.0) / (32.0 * std::pow(M_PI, 2.0) * physicalConstants::freeSpaceImpedance);
			auto freq_val = const_term * (std::pow(std::abs(L_pair.second + physicalConstants::freeSpaceImpedance * N_pair.first), 2.0) + std::pow(std::abs(L_pair.first - physicalConstants::freeSpaceImpedance * N_pair.second), 2.0));
			pot_rad_[angpair][frequencies[f]] = freq_val;
		}
	}
}

}
