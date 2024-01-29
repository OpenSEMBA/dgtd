#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;


void func_exp_real_part_2D(const Vector& x, Vector& v, const double freq, const Rho angle)
{
	v[0] = cos(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
	v[1] = cos(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
}

void func_exp_imag_part_2D(const Vector& x, Vector& v, const double freq, const Rho angle)
{
	v[0] = sin(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
	v[1] = sin(2.0 * M_PI * freq * sqrt(std::pow(x[0], 2.0) + std::pow(x[1], 2.0)) * cos(angle));
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

GridFunction assembleHigherVDimL2GridFunction(const GridFunction& xField, const GridFunction& yField)
{
	auto l2fes_vdim{
		FiniteElementSpace(
			xField.FESpace()->GetMesh(),
			dynamic_cast<const L2_FECollection*>(xField.FESpace()->FEColl()),
			2)
	};

	GridFunction res(&l2fes_vdim);
	res.SetVector(xField, 0);
	res.SetVector(yField, xField.Size());

	return res;
}

std::pair<Nedelec_XY, H1_Z> convertL2toNewFES2D(const GridFunction& xField, const GridFunction& yField, const GridFunction& zField)
{

	if (xField.FESpace() != yField.FESpace() || xField.FESpace() != zField.FESpace()) {
		throw std::exception("GridFunctions xField, yField and/or zField do not have the same FiniteElementSpace.");
	}
	auto l2_gf_vdim{ &assembleHigherVDimL2GridFunction(xField, yField) };
	VectorGridFunctionCoefficient dg_vgfc(l2_gf_vdim);
	auto ndfes{	FiniteElementSpace(	xField.FESpace()->GetMesh(), &ND_FECollection(xField.FESpace()->FEColl()->GetOrder(), xField.FESpace()->GetMesh()->Dimension())) };
	Nedelec_XY xyFieldND(&ndfes);
	xyFieldND.ProjectCoefficient(dg_vgfc);

	auto h1fes{ FiniteElementSpace( zField.FESpace()->GetMesh(), &H1_FECollection(zField.FESpace()->FEColl()->GetOrder(), zField.FESpace()->GetMesh()->Dimension())) };
	H1_Z zFieldH1(&h1fes);
	zFieldH1.ProjectGridFunction(zField);

	return std::make_pair(xyFieldND, zFieldH1);

}


std::unique_ptr<LinearForm> assembleLinearForm(VectorFunctionCoefficient& vfc, FiniteElementSpace& fes) 
{
	auto res{ std::make_unique<LinearForm>(&fes) };
	switch (vfc.GetVDim()) {
	case 1:
		res->AddBdrFaceIntegrator(new BoundaryTangentialLFIntegrator(vfc), getNTFFMarker(fes.GetMesh()->bdr_attributes.Max()));
		break;
	case 2:
		res->AddBdrFaceIntegrator(new VectorFEBoundaryTangentLFIntegrator(vfc), getNTFFMarker(fes.GetMesh()->bdr_attributes.Max()));
		break;
	default:
		throw std::runtime_error("Wrong VDIM for VectorFunctionCoefficient.");
	}
	res->Assemble();
	return res;
}


VectorFunctionCoefficient buildVFC_2D(const double freq, const Rho& angle, bool isReal)
{

	std::function<void(const Vector&, Vector&)> f = 0;
	switch (isReal) {
		case true:
			f = std::bind(&func_exp_real_part_2D, std::placeholders::_1, std::placeholders::_2, freq, angle);
			break;
		case false:
			f = std::bind(&func_exp_imag_part_2D, std::placeholders::_1, std::placeholders::_2, freq, angle);
			break;
	}
	VectorFunctionCoefficient res(2, f);
	return res;

}

void RCSManager::performRCS2DCalculations(GridFunction& Ax, GridFunction& Ay, GridFunction& Bz, const double frequency, const std::pair<Rho,Phi>& angles)
{
	auto pair = convertL2toNewFES2D(Ax, Ay, Bz);

	auto vfc_real_nd{ buildVFC_2D(frequency, angles.first, true) };
	auto lf_real_nd{ assembleLinearForm(vfc_real_nd, *pair.first.FESpace()) };
	//eval with gf, gives a double, add them all and map them to the freq and angles

	auto vfc_imag_nd{ buildVFC_2D(frequency, angles.first, false) };
	auto lf_imag_nd{ assembleLinearForm(vfc_imag_nd, *pair.first.FESpace()) };

	auto vfc_real_h1{ buildVFC_2D(frequency, angles.first, true) };
	auto lf_real_h1{ assembleLinearForm(vfc_real_h1, *pair.second.FESpace()) };

	auto vfc_imag_h1{ buildVFC_2D(frequency, angles.first, false) };
	auto lf_imag_h1{ assembleLinearForm(vfc_imag_h1, *pair.second.FESpace()) };

}


RCSManager::RCSManager(const std::string& path, const std::vector<double>& frequency, const SphericalAngles& angle) 
{
	basePath_ = path;

	std::unique_ptr<GridFunction> Ex, Ey, Ez, Hx, Hy, Hz;
	double time;
	Mesh mesh{ Mesh::LoadFromFile(basePath_ + "/mesh", 1, 0) };

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
			time = getTime(dir_entry.path().generic_string() + "/time.txt");

			for (const auto& f : frequency) {
				for (const auto& angpair : angle) {
					switch (mesh.SpaceDimension()) {
					case 2:
						performRCS2DCalculations(*Ex, *Ey, *Hz, f, angpair);
						performRCS2DCalculations(*Hx, *Hy, *Ez, f, angpair);
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
	
}



}