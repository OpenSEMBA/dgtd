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

void RCSManager::performRCS2DCalculations(const std::string& path) 
{
	auto pairTE = convertL2toNewFES2D(getGridFunction(path, E, X), getGridFunction(path, E, Y), getGridFunction(path, H, Z));
	auto pairTM = convertL2toNewFES2D(getGridFunction(path, H, X), getGridFunction(path, H, Y), getGridFunction(path, E, Z));

	VectorFunctionCoefficient vfc_te_real_nd(pairTE.first.VectorDim(), func_exp_real_part);
	auto lf_real_nd{ assembleLinearForm(vfc_te_real_nd, *pairTE.first.FESpace()) };

	VectorFunctionCoefficient vfc_te_imag_nd(pairTE.first.VectorDim(), func_exp_imag_part);
	auto lf_imag_nd{ assembleLinearForm(vfc_te_imag_nd, *pairTE.first.FESpace()) };

	VectorFunctionCoefficient vfc_te_real_h1(pairTE.second.VectorDim(), func_exp_real_part);
	auto lf_real_h1{ assembleLinearForm(vfc_te_real_h1, *pairTE.second.FESpace()) };

	VectorFunctionCoefficient vfc_te_imag_h1(pairTE.second.VectorDim(), func_exp_imag_part);
	auto lf_imag_h1{ assembleLinearForm(vfc_te_imag_h1, *pairTE.second.FESpace()) };

	//Same for TM. Duplicating code, gotta fix. Can probably economise down to single assembly with double coeff, or viceversa.
}

RCSManager::RCSManager(const std::string& path, const std::string& probe_name, double f) 
{
	basePath_ = path;
	frequency_ = f;
	//frequency_ = somethingSomethingFrequency();
	m_ = Mesh::LoadFromFile(path + "/" + probe_name + "_000000/mesh", 1, 0);
	for (auto const& dir_entry : std::filesystem::directory_iterator{ basePath_ }) {
		auto iterPath{ dir_entry.path().generic_string() };
		switch (getGridFunction(path, E, X).FESpace()->GetMesh()->Dimension()) {
		case 2:
			performRCS2DCalculations(iterPath);
			break;
		case 3:
			break;
		default:
			throw std::runtime_error("The RCSManager does not support dimensions other than 2 or 3.");
		}
	}
	
}

GridFunction RCSManager::getGridFunction(const std::string& path, const FieldType& f, const Direction& d)
{
	auto filePath{ getGridFunctionPathForType(path, f, d) };
	std::ifstream in(filePath);
	if (!in) {
		throw std::runtime_error("File could not be opened in getGridFunctionPathForType.");
	}
	GridFunction gf(&m_, in);
	return gf;
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