#include "OpensembaAdapter.h"

#include "parsers/json/Parser.h"
#include "core/ProblemDescription.h"
#include "core/util/OptionsBase.h"
#include "core/physicalModel/physicalModels.h"

using semba::util::OptionsBase;
using semba::PMGroup;
using json = nlohmann::json;


namespace maxwell {

FluxType strToFluxType(const std::string& label)
{
	if (label == "centered") {
		return FluxType::Centered;
	}
	else if (label == "upwind") {
		return FluxType::Upwind;
	}
	else {
		throw std::runtime_error("Invalid flux type label.");
	}
}

FieldType strToFieldType(const std::string& label)
{
	if (label == "E") {
		return FieldType::E;
	}
	else if (label == "H") {
		return FieldType::H;
	}
	else {
		throw std::runtime_error("Invalid field type label.");
	}
}

template <class T>
std::vector<T> readVector(const json& j)
{
	std::vector<T> r;
	for (const auto& v : j) {
		r.push_back(v.get<T>());
	}
	return r;
}

mfem::Vector toMFEMVector(const std::vector<double>& v)
{
	mfem::Vector r(v.size());
	std::copy(v.begin(), v.end(), r.begin());
	return r;
}

OpensembaAdapter::OpensembaAdapter(const std::string& fn) :
	filename_{fn}
{
	std::ifstream stream(filename_);
	if (!stream.is_open()) {
		throw std::logic_error("Can not open file: " + filename_);
	}

	try {
		stream >> json_;
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}
}

std::unique_ptr<Function> readFunction(const json& j)
{
	auto type{ j["type"].get<std::string>() };
	if (type == "sinusoidalMode") {
		return std::make_unique<SinusoidalMode>(
			readVector<std::size_t>(j["modes"])
		);
	}
	else {
		throw std::runtime_error("Unsupported function.");
	}
}

InitialField readInitialFieldSource(const json& j)
{
	return {
		*readFunction(j["function"]),
		strToFieldType(j["fieldType"].get<std::string>()),
		toMFEMVector(readVector<double>(j["polarization"])),
		toMFEMVector(readVector<double>(j["position"])),
	};
}

Sources readSources(const json& j)
{
	auto srcs{ j.find("sources")};
	Sources r;
	for (const auto& s : *srcs) {
		auto type{ s["sourceType"].get<std::string>() };
		if (type == "initialField") {
			r.add(readInitialFieldSource(s));
		}
		else {
			throw std::runtime_error(
				"Unsupported source type."
			);
		}
	}
	return r;
}

ExporterProbe readExporterProbe(const json& j)
{
	return {
		j["name"].get<std::string>(),
		j["visSteps"].get<int>()
	};
}

Probes readProbes(const json& j)
{
	auto ps{ j.find("probes") };
	Probes r;
	for (const auto& p : *ps) {
		auto type{ p["type"].get<std::string>() };
		if (type == "exporter") {
			r.exporterProbes.push_back(readExporterProbe(p));
		}
		else {
			throw std::runtime_error(
				"Unsupported probe type."
			);
		}
	}
	return r;
}

mfem::Mesh OpensembaAdapter::readMesh(const json& j) const
{
	if (j["mesh"]["type"] != "gmsh") {
		throw std::runtime_error(
			"Invalid mesh type. Only gmsh files are supported"
		);
	}

	semba::util::ProjectFile fn{ filename_ };
	auto caseName{
		semba::util::ProjectFile::removeExtension(
			semba::util::ProjectFile::removeExtension(
				fn.getBasename()
			)
		)
	};
	auto meshFilename{fn.getFolder() + caseName + ".msh" };
	return mfem::Mesh{ meshFilename.c_str() };
}

Model OpensembaAdapter::readModel(const json& j) const
{
	auto modelJSON{ j.find("model") };
	if (modelJSON == j.end()) {
		throw std::runtime_error("Can not find \"model\" label.");
	}

	AttributeToMaterial attToMat;
	AttributeToBoundary attToBdr;
	AttributeToBoundary attToInteriorConditions;
	for (const auto& m : semba::parsers::JSON::readMaterials(*modelJSON)) {
		Attribute att( m->getId().toInt() );
		if (m->is<semba::physicalModel::Vacuum>()) {
			attToMat.emplace(att, buildVacuumMaterial());
		} 
		else if (m->is<semba::physicalModel::PEC>()) {
			attToBdr.emplace(att, BdrCond::PEC);
		}
		else {
			throw std::runtime_error("Unsupported material.");
		}
	}

	return { 
		readMesh(*modelJSON), 
		attToMat, 
		attToBdr, 
		attToInteriorConditions 
	};
}

Problem OpensembaAdapter::readProblem() const
{
	try {
		return {
			readModel(json_) ,
			readProbes(json_),
			readSources(json_)
		};
	}
	catch (json::exception& e) {
		std::cout << e.what() << '\n';
	}
	
}

EvolutionOptions readEvolutionOptions(const json& j)
{
	EvolutionOptions r;
	OptionsBase::setIfExists<int>(j, r.order, "order");
	OptionsBase::setIfExists<bool>(j, r.spectral,"spectral");
	
	OptionsBase::setIfExistsUsingLabelConversion<FluxType>(
		j, r.fluxType, "fluxType", strToFluxType);

	return r;
}

SolverOptions OpensembaAdapter::readSolverOptions() const
{
	auto& analysis{json_.find("analysis")};
	if (analysis == json_.end()) {
		throw std::logic_error("File does not contain \"analysis\".");
	}

	if (analysis->at("solver").get<std::string>() != "opensemba/dgtd") {
		throw std::logic_error("Invalid solver input.");
	}

	if (analysis->at("units").get<std::string>() != "natural") {
		throw std::logic_error("Invalid units.");
	}

	auto& sO{ analysis->at("solverOptions") };

	SolverOptions r;
    OptionsBase::setIfExists<double>(sO, r.finalTime, "finalTime");
	OptionsBase::setIfExists<double>(sO, r.timeStep,  "timeStep");
	OptionsBase::setIfExists<double>(sO, r.cfl,       "cfl");
	
	auto& evolution{ sO.find("evolution")};
	if (evolution != sO.end()) {
		r.evolution = readEvolutionOptions(*evolution);
	}
	
	return r;
};


}