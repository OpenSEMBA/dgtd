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

void RCSManager::update(FieldsToTime& ftt)
{
	int rank{ 0 };

	for (const auto& itEntry : std::filesystem::directory_iterator(basePath_)) {
		auto subPath{ itEntry.path() };
		auto time{ getTime(subPath.generic_string() + "/time.txt") };
		for (auto& f : { E, H }) {
			for (auto& d : { X, Y, Z }) {
				calculateRCS(fieldsRCS_.at(getGridFunctionString(f, d)), 
					getGridFunction(subPath.string(), f, d), 
					getTime(subPath.generic_string() + "/time.txt")
				);
			}
		}
		rank++;
	}
	for (auto& f : { E, H }) {
		for (auto& d : { X, Y, Z }) {
			auto vec{ fieldsRCS_.at(getGridFunctionString(f, d)) };
			for (auto i{ 0 }; i < vec.size(); ++i) {
				vec[i] /= (double) rank;
			}
		}
	}
	//Iterate through directory X
	//Read GF data and assemble ftt X
	//Pass ftt to calculate X
	//Rinse and repeat while there's folders to iterate through
	//Keep counting how many folder we're going through, at the end, divide by that.
	
}
 
void RCSManager::calculateRCS(CompVec& rcs, const GridFunction& ingf, const Time time)
{
	const std::complex<double> constPart{ 0.0, -(double)(2.0 * M_PI * frequency_) };
	std::complex<double> res{ 0.0, 0.0 };
	for (auto i{ 0 }; i < ingf.Size(); ++i) {
		rcs[i] = ingf[i] * exp(constPart * time);
	}
	//Do the math
	//rcs.atspecificthing += mathmath(gf);

}

void RCSManager::initFieldsRCS(const std::string& path)
{
	for (auto f : { E, H }) {
		for (auto d : { X, Y, Z }) {
			fieldsRCS_.emplace(getGridFunctionString(f, d), getGridFunction(path, f, d));
		}
	}
}


}