#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;

double speed_of_light{ 299792458.0 };

RCSData::RCSData(double val, double f, SphericalAngles angles, double t) :
	RCSvalue{ val },
	frequency{ f },
	angles{ angles },
	time{ t }
{}

double func_exp_real_part_2D(const Vector& x, const double freq, const Rho angle)
{
	return cos(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
}

double func_exp_imag_part_2D(const Vector& x, const double freq, const Rho angle)
{
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


FunctionCoefficient buildFC_2D(const double freq, const Rho& angle, bool isReal)
{
	std::function<double(const Vector&)> f = 0;
	switch (isReal) {
		case true:
			f = std::bind(&func_exp_real_part_2D, std::placeholders::_1, freq, angle);
			break;
		case false:
			f = std::bind(&func_exp_imag_part_2D, std::placeholders::_1, freq, angle);
			break;
	}
	FunctionCoefficient res(f);
	return res;

}

double RCSManager::performRCS2DCalculations(GridFunction& Ax, GridFunction& Ay, GridFunction& Az, const double frequency, const SphericalAngles& angles)
{
	auto fc_real_l2{ buildFC_2D(frequency, angles.first, true) };
	auto lf_real_l2_x{ assembleLinearForm(fc_real_l2, *Ax.FESpace(), X) };
	auto lf_real_l2_y{ assembleLinearForm(fc_real_l2, *Ax.FESpace(), Y) };
	auto lf_real_l2_z{ assembleLinearForm(fc_real_l2, *Ax.FESpace(), Z) };

	auto fc_imag_l2{ buildFC_2D(frequency, angles.first, false) };
	auto lf_imag_l2_x{ assembleLinearForm(fc_imag_l2, *Ax.FESpace(), X) };
	auto lf_imag_l2_y{ assembleLinearForm(fc_imag_l2, *Ax.FESpace(), Y) };
	auto lf_imag_l2_z{ assembleLinearForm(fc_imag_l2, *Ax.FESpace(), Z) };

	//Manual extension of n x A = (nyAz - nzAy) dir_x + (nzAx - nxAz) dir_y + (nxAy - nyAx) dir_z

	auto real_factor = mfem::InnerProduct(*lf_real_l2_y.get(), Az) - mfem::InnerProduct(*lf_real_l2_z.get(), Ay)
					 + mfem::InnerProduct(*lf_real_l2_z.get(), Ax) - mfem::InnerProduct(*lf_real_l2_x.get(), Az)
					 + mfem::InnerProduct(*lf_real_l2_x.get(), Ay) - mfem::InnerProduct(*lf_real_l2_y.get(), Ax);

	auto imag_factor = mfem::InnerProduct(*lf_imag_l2_y.get(), Az) - mfem::InnerProduct(*lf_imag_l2_z.get(), Ay)
					 + mfem::InnerProduct(*lf_imag_l2_z.get(), Ax) - mfem::InnerProduct(*lf_imag_l2_x.get(), Az)
					 + mfem::InnerProduct(*lf_imag_l2_x.get(), Ay) - mfem::InnerProduct(*lf_imag_l2_y.get(), Ax);

	return sqrt(std::pow(real_factor, 2.0) + std::pow(imag_factor, 2.0));
}

void RCSManager::fillPostDataMaps(const std::vector<double>& frequencies, const SphericalAnglesVector& angleVec)
{
	for (const auto& angpair : angleVec) {
		Freq2Value inner;
		for (const auto& f : frequencies) {
			inner.emplace(f, 0.0);
		}
		postdata_.emplace(angpair, inner);
	}
}

RCSManager::RCSManager(const std::string& path, const std::vector<double>& frequencies, const SphericalAnglesVector& angleVec)
{
	basePath_ = path;

	auto auxiliarySpaces = false;

	std::unique_ptr<GridFunction> Ex, Ey, Ez, Hx, Hy, Hz;
	double time;
	int time_steps{ 0 };
	Mesh mesh{ Mesh::LoadFromFile(basePath_ + "/mesh", 1, 0) };

	fillPostDataMaps(frequencies, angleVec);

	for (auto const& dir_entry : std::filesystem::directory_iterator(basePath_)) {
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) != "mesh") {
			std::ifstream inEx(dir_entry.path().generic_string() + "/Ex.gf");
			Ex = std::make_unique<GridFunction>(&mesh, inEx);
			std::ifstream inEy(dir_entry.path().generic_string() + "/Ey.gf");
			Ey = std::make_unique<GridFunction>(&mesh, inEy);
			std::ifstream inEz(dir_entry.path().generic_string() + "/Ez.gf");
			Ez = std::make_unique<GridFunction>(&mesh, inEz);
			std::ifstream inHx(dir_entry.path().generic_string() + "/Hx.gf");
			Hx = std::make_unique<GridFunction>(&mesh, inHx);
			std::ifstream inHy(dir_entry.path().generic_string() + "/Hy.gf");
			Hy = std::make_unique<GridFunction>(&mesh, inHy);
			std::ifstream inHz(dir_entry.path().generic_string() + "/Hz.gf");
			Hz = std::make_unique<GridFunction>(&mesh, inHz);
			time = getTime(dir_entry.path().generic_string() + "/time.txt") / speed_of_light;
			time_steps++;

			std::unique_ptr<RCSData> data;
			for (const auto& f : frequencies) {
				for (const auto& angpair : angleVec) {
					switch (mesh.SpaceDimension()) {
					case 2:					
						data = std::make_unique<RCSData>(performRCS2DCalculations(*Ex, *Ey, *Ez, f, angpair) + performRCS2DCalculations(*Hx, *Hy, *Hz, f, angpair), f, angpair, time);
						postdata_[angpair][f] += data->RCSvalue;
						break;
					case 3:
						break;
					default:
						throw std::runtime_error("RCSManager cannot be applied on dimensions other than 2 and 3.");
					}
				}
			}
		}
	}

	for (const auto& angpair : angleVec) {
		std::ofstream myfile;
		myfile.open("RCSData_" + std::to_string(angpair.first) + "_" + std::to_string(angpair.second) + "_dgtd.dat");	
		myfile << "Angle Rho " << "Angle Phi " << "Frequency (Hz) " << "10*log(RCSData/Lambda)\n";
		for (const auto& f : frequencies) {
			postdata_[angpair][f] /= (time_steps * 2.0); //( 2.0 = Correction value for 2D)
			myfile << angpair.first << " " << angpair.second << " " << f << " " << 10.0*log(postdata_[angpair][f]*f) << "\n";			
		}
		myfile.close();
	}

	
}


}