#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;
using Nedelec_XY = GridFunction;
using H1_Z = GridFunction;

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
	auto h1fes{ FiniteElementSpace( zField.FESpace()->GetMesh(), &H1_FECollection(zField.FESpace()->FEColl()->GetOrder(), zField.FESpace()->GetMesh()->Dimension())) };
	
	Nedelec_XY xyFieldND(&ndfes);
	xyFieldND.ProjectCoefficient(dg_vgfc);
	H1_Z zFieldH1(&h1fes);
	zFieldH1.ProjectGridFunction(zField);

	return std::make_pair(xyFieldND, zFieldH1);

}

Mesh getRCSMesh(const std::string& path)
{
	std::ifstream in(path + "/mesh");
	return Mesh(in);
}

std::string getGridFunctionString(const FieldType& f, const Direction& d)
{
	switch (f) {
	case E:
		switch (d) {
		case X:
			return "Ex";
		case Y:
			return "Ey";
		case Z:
			return "Ez";
		}
	case H:
		switch (d) {
		case X:
			return "Hx";
		case Y:
			return "Hy";
		case Z:
			return "Hz";
		}
	}
}

std::string getGridFunctionPathForType(const std::string& path, const FieldType& f, const Direction& d)
{
	return path + "/" + getGridFunctionString(f, d) + ".gf";
}

void func_exp_real_part(const Vector& pos, double time, Vector& v)
{
	//real part of exponential things
}

void func_exp_imag_part(const Vector& pos, double time, Vector& v)
{
	//real part of exponential things
}

Array<int> getNTFFMarker(const int att_size)
{
	Array<int> res(att_size);
	res = 0;
	res[static_cast<int>(BdrCond::NearToFarField) - 1] = 1;
	return res;
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

void RCSManager::performRCS2DCalculations(GridFunction& Ax, GridFunction& Ay, GridFunction& Bz) 
{
	auto pairTE = convertL2toNewFES2D(Ax, Ay, Bz);

	VectorFunctionCoefficient vfc_te_real_nd(pairTE.first.VectorDim(), func_exp_real_part);
	auto lf_real_nd{ assembleLinearForm(vfc_te_real_nd, *pairTE.first.FESpace()) };

	VectorFunctionCoefficient vfc_te_imag_nd(pairTE.first.VectorDim(), func_exp_imag_part);
	auto lf_imag_nd{ assembleLinearForm(vfc_te_imag_nd, *pairTE.first.FESpace()) };

	VectorFunctionCoefficient vfc_te_real_h1(pairTE.second.VectorDim(), func_exp_real_part);
	auto lf_real_h1{ assembleLinearForm(vfc_te_real_h1, *pairTE.second.FESpace()) };

	VectorFunctionCoefficient vfc_te_imag_h1(pairTE.second.VectorDim(), func_exp_imag_part);
	auto lf_imag_h1{ assembleLinearForm(vfc_te_imag_h1, *pairTE.second.FESpace()) };



}

RCSManager::RCSManager(const std::string& path, const std::string& probe_name, double f) 
{
	basePath_ = path;
	frequency_ = f;

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

			switch (mesh.SpaceDimension()) {
			case 2:
				performRCS2DCalculations(*Ex, *Ey, *Hz);
				performRCS2DCalculations(*Hx, *Hy, *Ez);
				break;
			case 3:
				break;
			default:
				throw std::runtime_error("RCSManager cannot be applied on dimensions other than 2 and 3.");
			}

		}
	}
	
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



}