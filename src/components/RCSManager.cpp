#include <components/RCSManager.h>

namespace maxwell {

using namespace mfem;

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
			return "Ex.gf";
		case Y:
			return "Ey.gf";
		case Z:
			return "Ez.gf";
		}
	case H:
		switch (d) {
		case X:
			return "Hx.gf";
		case Y:
			return "Hy.gf";
		case Z:
			return "Hz.gf";
		}
	}
}

std::string getGridFunctionPathForType(const std::string& path, const FieldType& f, const Direction& d)
{
	return path + "/" + getGridFunctionString(f, d);
}

RCSManager::RCSManager(const std::string& path, const NearToFarFieldProbe& p) 
{
	basePath_ = path;
	m_ = std::make_unique<Mesh>(getRCSMesh(basePath_));
	const std::string initFolder{ path + "/" + p.name + "_000000" };
	initFieldsRCS(initFolder);
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

void RCSManager::update(FieldsToTime& ftt)
{
	
}

void RCSManager::initFieldsRCS(const std::string& path)
{
	for (auto f : { E,H }) {
		for (auto d : { X, Y, Z }) {
			fieldsRCS_.emplace(getGridFunctionString(f, d), getGridFunction(path, f, d));
		}
	}
}

void RCSManager::calculateRCS(FieldsToTime& ftt)
{

}


}