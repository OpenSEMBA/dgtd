#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;

double speed_of_light{ 299792458.0 };
double mu_0 = 4.0e-7 * M_PI;
double e_0 = 8.8541878128e-12;
double eta_0{ sqrt(mu_0 / e_0) };

double func_exp_real_part_2D(const Vector& x, const double freq, const Phi phi)
{
	//angulo viene dado por x[0], x[1] y 0.0, 0.0. No es el angulo donde observo, es el angulo que forma el punto y el angulo de observacion en un sistema centrado en el punto.
	auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
	return cos(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
}

double func_exp_imag_part_2D(const Vector& x, const double freq, const Phi phi)
{
	auto angle = acos((x[0] * cos(phi) + x[1] * sin(phi)) / sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)));
	return sin(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
}

Array<int> getNTFFMarker(const int att_size)
{
	Array<int> res(att_size);
	res = 0;
	res[static_cast<int>(BdrCond::NearToFarField) - 1] = 1;
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

GridFunction parseGridFunction(Mesh& mesh, const std::string& path)
{
	std::ifstream in(path);
	GridFunction res(&mesh, in);
	return res;
}

DFTFreqFieldsComplex RCSManager::assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& field)
{
	DFTFreqFieldsComplex res(frequencies.size());
	size_t dofs{ 0 };
	auto time_step_counter{ 0 };
	double time_step;
	std::vector<double> two_times_for_dt(2);
	size_t time_counter{ 0 };

	for (auto const& dir_entry : std::filesystem::directory_iterator(base_path_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			two_times_for_dt[time_counter] = getTime(dir_entry.path().generic_string() + "/time.txt");
			if (time_counter == 1) {
				break;
			}
			time_counter++;
		}
	}

	time_step = std::abs(two_times_for_dt[1] - two_times_for_dt[0]);

	for (auto const& dir_entry : std::filesystem::directory_iterator(base_path_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			auto gf{ parseGridFunction(mesh, dir_entry.path().generic_string() + field) };
			auto time = getTime(dir_entry.path().generic_string() + "/time.txt");
			dofs = gf.Size();
			time_step_counter++;
			for (int f{ 0 }; f < frequencies.size(); f++) {
				auto time_const{ time_step * std::exp(std::complex<double>(0.0, -2.0 * M_PI * frequencies[f] * time)) };
				res[f].resize(dofs);
				for (int i{ 0 }; i < dofs; i++) {
					res[f][i] += std::complex<double>(gf[i], 0.0) * time_const;
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

void splitCompIntoDoubles(const std::vector<std::complex<double>>& comp, Vector& real, Vector& imag)
{
	for (int i{ 0 }; i < comp.size(); i++) {
		real[i] = comp[i].real();
		imag[i] = comp[i].imag();
	}
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

std::pair<std::complex<double>, std::complex<double>> RCSManager::performRCS2DCalculations(ComplexVector& FAx, ComplexVector& FAy, ComplexVector& FAz, const double frequency, const SphericalAngles& angles, bool isElectric)
{
	auto fc_real_l2{ buildFC_2D(frequency, angles.first, true) };
	auto fc_imag_l2{ buildFC_2D(frequency, angles.first, false) };

	std::unique_ptr<LinearForm> lf_real_l2_x, lf_real_l2_y, lf_real_l2_z, lf_imag_l2_x, lf_imag_l2_y, lf_imag_l2_z;
	switch (isElectric) {
	case true:
		lf_real_l2_x = assembleLinearForm(fc_real_l2, *gfex_->FESpace(), X);
		lf_real_l2_y = assembleLinearForm(fc_real_l2, *gfey_->FESpace(), Y);
		lf_real_l2_z = assembleLinearForm(fc_real_l2, *gfez_->FESpace(), Z);
		lf_imag_l2_x = assembleLinearForm(fc_imag_l2, *gfex_->FESpace(), X);
		lf_imag_l2_y = assembleLinearForm(fc_imag_l2, *gfey_->FESpace(), Y);
		lf_imag_l2_z = assembleLinearForm(fc_imag_l2, *gfez_->FESpace(), Z);
		break;
	case false:
		lf_real_l2_x = assembleLinearForm(fc_real_l2, *gfhx_->FESpace(), X);
		lf_real_l2_y = assembleLinearForm(fc_real_l2, *gfhy_->FESpace(), Y);
		lf_real_l2_z = assembleLinearForm(fc_real_l2, *gfhz_->FESpace(), Z);
		lf_imag_l2_x = assembleLinearForm(fc_imag_l2, *gfhx_->FESpace(), X);
		lf_imag_l2_y = assembleLinearForm(fc_imag_l2, *gfhy_->FESpace(), Y);
		lf_imag_l2_z = assembleLinearForm(fc_imag_l2, *gfhz_->FESpace(), Z);
		break;
	}

	ComplexVector lf_x(lf_real_l2_x->Size()), lf_y(lf_real_l2_y->Size()), lf_z(lf_real_l2_z->Size());

	for (int i{ 0 }; i < lf_x.size(); i++) {
		lf_x[i] = std::complex<double>(lf_real_l2_x->Elem(i), lf_imag_l2_x->Elem(i));
		lf_y[i] = std::complex<double>(lf_real_l2_y->Elem(i), lf_imag_l2_y->Elem(i));
		lf_z[i] = std::complex<double>(lf_real_l2_z->Elem(i), lf_imag_l2_z->Elem(i));
	}

	auto DCx{ complexInnerProduct(lf_y, FAz) - complexInnerProduct(lf_z, FAy) };
	auto DCy{ complexInnerProduct(lf_z, FAx) - complexInnerProduct(lf_x, FAz) };
	auto DCz{ complexInnerProduct(lf_x, FAy) - complexInnerProduct(lf_y, FAx) };

	auto phi_value = -DCx * sin(angles.first)                      + DCy * cos(angles.first);
	auto rho_value =  DCx * cos(angles.second) * cos(angles.first) + DCy * cos(angles.second) * sin(angles.first) - DCz * sin(angles.second);

	return std::pair<std::complex<double>, std::complex<double>>(phi_value, rho_value);
}

RCSManager::RCSManager(const std::string& path, const std::vector<double>& frequencies, const std::vector<SphericalAngles>& angle_vec)
{
	base_path_ = path;
	
	Mesh mesh{ Mesh::LoadFromFile(base_path_ + "/mesh", 1, 0) };

	for (auto const& dir_entry : std::filesystem::directory_iterator(base_path_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			std::ifstream inex(dir_entry.path().generic_string() + "/Ex.gf");
			gfex_ = std::make_unique<GridFunction>(&mesh, inex);
			std::ifstream iney(dir_entry.path().generic_string() + "/Ey.gf");
			gfey_ = std::make_unique<GridFunction>(&mesh, iney);
			std::ifstream inez(dir_entry.path().generic_string() + "/Ez.gf");
			gfez_ = std::make_unique<GridFunction>(&mesh, inez);
			std::ifstream inhx(dir_entry.path().generic_string() + "/Hx.gf");
			gfhx_ = std::make_unique<GridFunction>(&mesh, inhx);
			std::ifstream inhy(dir_entry.path().generic_string() + "/Hy.gf");
			gfhy_ = std::make_unique<GridFunction>(&mesh, inhy);
			std::ifstream inhz(dir_entry.path().generic_string() + "/Hz.gf");
			gfhz_ = std::make_unique<GridFunction>(&mesh, inhz);
			break;
		}
	}

	fillPostDataMaps(frequencies, angle_vec);
	
	DFTFreqFieldsComplex FEx{ assembleFreqFields(mesh, frequencies, "/Ex.gf") };
	DFTFreqFieldsComplex FEy{ assembleFreqFields(mesh, frequencies, "/Ey.gf") };
	DFTFreqFieldsComplex FEz{ assembleFreqFields(mesh, frequencies, "/Ez.gf") };
 	DFTFreqFieldsComplex FHx{ assembleFreqFields(mesh, frequencies, "/Hx.gf") };
	DFTFreqFieldsComplex FHy{ assembleFreqFields(mesh, frequencies, "/Hy.gf") };
	DFTFreqFieldsComplex FHz{ assembleFreqFields(mesh, frequencies, "/Hz.gf") };

	std::unique_ptr<RCSData> data;
	double freqdata;
	double const_term;
	std::pair<std::complex<double>, std::complex<double>> N_pair, L_pair;
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			switch (mesh.SpaceDimension()) {
			case 2:	
				N_pair = performRCS2DCalculations(FEx[f], FEy[f], FEz[f], frequencies[f], angpair, true);
				L_pair = performRCS2DCalculations(FHx[f], FHy[f], FHz[f], frequencies[f], angpair, false);
				const_term = std::pow((2.0 * M_PI * frequencies[f]), 2.0) / (8.0 * M_PI * eta_0);
				freqdata = const_term * (std::pow(std::norm(L_pair.first + eta_0 * N_pair.second), 2.0) + std::pow(std::norm(L_pair.second - eta_0 * N_pair.first), 2.0));
				data = std::make_unique<RCSData>(freqdata, frequencies[f], angpair);
				postdata_[angpair][frequencies[f]] += data->RCSvalue;
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
		myfile.open("RCSData_" + std::to_string(angpair.first) + "_" + std::to_string(angpair.second) + "_dgtd.dat");
		myfile << "Angle Rho " << "Angle Phi " << "Frequency (Hz) " << "10*log(RCSData/Lambda)\n";
		for (const auto& f : frequencies) {
			myfile << angpair.first << " " << angpair.second << " " << f << " " << 10.0 * log(postdata_[angpair][f] * f) << "\n";
		}
		myfile.close();
	}

}


}