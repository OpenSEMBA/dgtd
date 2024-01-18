#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;
using NedelecXY = GridFunction;
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

std::pair<NedelecXY, H1_Z> convertL2toNewFES2D(const GridFunction& xField, const GridFunction& yField, const GridFunction& zField)
{

	if (xField.FESpace() != yField.FESpace() || xField.FESpace() != zField.FESpace()) {
		throw std::exception("GridFunctions xField, yField and/or zField do not have the same FiniteElementSpace.");
	}
	auto l2_gf_vdim{ &assembleHigherVDimL2GridFunction(xField, yField) };
	VectorGridFunctionCoefficient dg_vgfc(l2_gf_vdim);

	auto ndfes{	FiniteElementSpace(	xField.FESpace()->GetMesh(), &ND_FECollection(xField.FESpace()->FEColl()->GetOrder(), xField.FESpace()->GetMesh()->Dimension())) };
	auto h1fes{ FiniteElementSpace( zField.FESpace()->GetMesh(), &H1_FECollection(zField.FESpace()->FEColl()->GetOrder(), zField.FESpace()->GetMesh()->Dimension())) };
	
	NedelecXY xyFieldND(&ndfes);
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

RCSManager::RCSManager(const std::string& path, const NearToFarFieldProbe& p) 
{
	basePath_ = path;
	m_ = std::make_unique<Mesh>(getRCSMesh(basePath_));
	//frequency_ = somethingSomethingFrequency();
	const std::string initFolder{ path + "/" + p.name + "_000000" };
	
}

GridFunction RCSManager::getGridFunction(const std::string& path, const FieldType& f, const Direction& d)
{
	auto filepath{ getGridFunctionPathForType(path, f, d) };
	std::ifstream in(filepath);
	if (!in) {
		throw std::runtime_error("File could not be opened in readGridFunctionFromFile, verify path.");
	}
	GridFunction gf(m_.get(), in);
	return gf;
}

const double getTime(const std::string& timePath)
{
	std::ifstream timeFile(timePath);
	if (!timeFile) {
		throw std::runtime_error("File could not be opened in getTime for RCS, verify path.");
	}
	std::string timeString;
	std::getline(timeFile, timeString);
	return std::stod(timeString);
}



}