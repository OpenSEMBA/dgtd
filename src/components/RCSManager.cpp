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
	int rank{ 0 };

	for (auto itEntry = std::filesystem::recursive_directory_iterator(basePath_);
		itEntry != std::filesystem::recursive_directory_iterator();
		++itEntry) {



	}
	//Iterate through directory
	//Read GF data and assemble ftt
	//Pass ftt to calculate
	calculateRCS(ftt);
	//Rinse and repeat while there's folders to iterate through
	//Keep counting how many folder we're going through, at the end, divide by that.
	rank++;
}

void RCSManager::initFieldsRCS(const std::string& path)
{
	for (auto f : { E, H }) {
		for (auto d : { X, Y, Z }) {
			fieldsRCS_.emplace(getGridFunctionString(f, d), getGridFunction(path, f, d));
		}
	}
}

void RCSManager::calculateRCS(const GridFunction& rcs, GridFunction& gf)
{
	const std::complex<double> constPart{ 0.0, -(double)(2.0 * M_PI * frequency_) };
	//Do the math
	//rcs.atspecificthing += mathmath(gf);

}


}