#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;

double speed_of_light{ 299792458.0 };

RCSData::RCSData(double val, double f, SphericalAngles angles) :
	RCSvalue{ val },
	frequency{ f },
	angles{ angles }
{}

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

DFTFreqFieldsComp RCSManager::assembleFreqFields(Mesh& mesh, const std::vector<double>& frequencies, const std::string& field)
{
	DFTFreqFieldsComp res(frequencies.size());
	size_t dofs{ 0 };
	auto time_step_counter{ 0 };
	double time_step;
	std::vector<double> two_times_for_dt(2);
	size_t time_counter{ 0 };

	for (auto const& dir_entry : std::filesystem::directory_iterator(base_path_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			two_times_for_dt[time_counter] = getTime(dir_entry.path().generic_string() + "/time.txt") / speed_of_light;
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
			auto time = getTime(dir_entry.path().generic_string() + "/time.txt") / speed_of_light;
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

std::complex<double> complexInnerProduct(ComplexVector& first, ComplexVector& last)
{
	if (first.size() != last.size()) {
		throw std::runtime_error("Complex Vectors do not have the same sizes.");
	}
	std::complex<double> res(0.0, 0.0);
	for (int i{ 0 }; i < first.size(); i++) {
		res += first[i] * std::conj(last[i]);
	}
	return res;
}

std::pair<std::complex<double>, std::complex<double>> RCSManager::performRCS2DCalculations(ComplexVector& FAx, ComplexVector& FAy, ComplexVector& FAz, const double frequency, const SphericalAngles& angles)
{
	auto fc_real_l2{ buildFC_2D(frequency, angles.first, true) };
	auto lf_real_l2_x{ assembleLinearForm(fc_real_l2, *fes_, X) };
	auto lf_real_l2_y{ assembleLinearForm(fc_real_l2, *fes_, Y) };
	auto lf_real_l2_z{ assembleLinearForm(fc_real_l2, *fes_, Z) };

	auto fc_imag_l2{ buildFC_2D(frequency, angles.first, false) };
	auto lf_imag_l2_x{ assembleLinearForm(fc_imag_l2, *fes_, X) };
	auto lf_imag_l2_y{ assembleLinearForm(fc_imag_l2, *fes_, Y) };
	auto lf_imag_l2_z{ assembleLinearForm(fc_imag_l2, *fes_, Z) };

	ComplexVector lf_x(lf_real_l2_x->Size()), lf_y(lf_real_l2_y->Size()), lf_z(lf_real_l2_z->Size());

	for (int i{ 0 }; i < lf_x.size(); i++) {
		lf_x[i] = std::complex<double>(lf_real_l2_x.get()->Elem(i), lf_imag_l2_x.get()->Elem(i));
		lf_y[i] = std::complex<double>(lf_real_l2_y.get()->Elem(i), lf_imag_l2_y.get()->Elem(i));
		lf_z[i] = std::complex<double>(lf_real_l2_z.get()->Elem(i), lf_imag_l2_z.get()->Elem(i));
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
			std::ifstream in(dir_entry.path().generic_string() + "/Ex.gf");
			GridFunction gf(&mesh, in);
			fec_ = std::make_unique<DG_FECollection>(gf.FESpace()->FEColl()->GetOrder(), mesh.SpaceDimension(), 1);
			fes_ = std::make_unique<FiniteElementSpace>(&mesh, fec_.get());
			break;
		}
	}

	fillPostDataMaps(frequencies, angle_vec);
	
	DFTFreqFieldsComp FEx{ assembleFreqFields(mesh, frequencies, "/Ex.gf") };
	DFTFreqFieldsComp FEy{ assembleFreqFields(mesh, frequencies, "/Ey.gf") };
	DFTFreqFieldsComp FEz{ assembleFreqFields(mesh, frequencies, "/Ez.gf") };
 	DFTFreqFieldsComp FHx{ assembleFreqFields(mesh, frequencies, "/Hx.gf") };
	DFTFreqFieldsComp FHy{ assembleFreqFields(mesh, frequencies, "/Hy.gf") };
	DFTFreqFieldsComp FHz{ assembleFreqFields(mesh, frequencies, "/Hz.gf") };

	std::unique_ptr<RCSData> data;
	double freqdata;
	double const_term;
	std::pair<std::complex<double>, std::complex<double>> N_pair, L_pair;
	for (int f{ 0 }; f < frequencies.size(); f++) {
		for (const auto& angpair : angle_vec) {
			switch (mesh.SpaceDimension()) {
			case 2:	
				N_pair = performRCS2DCalculations(FEx[f], FEy[f], FEz[f], frequencies[f], angpair);
				L_pair = performRCS2DCalculations(FHx[f], FHy[f], FHz[f], frequencies[f], angpair);
				const_term = std::pow((2.0 * M_PI * frequencies[f]), 2.0) / (8.0 * M_PI);
				freqdata = const_term * (std::pow(std::norm(L_pair.first + N_pair.second), 2.0) + std::pow(std::norm(L_pair.second - N_pair.first), 2.0));
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