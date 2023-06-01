#include "OpensembaAdapter.h"

#include "parsers/json/Parser.h"
#include "core/ProblemDescription.h"
#include "core/util/OptionsBase.h"

using semba::util::OptionsBase;
using semba::PMGroup;
using json = nlohmann::json;

namespace maxwell {

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

//Sources readSources(const json& j)
//{
//	// TODO
//}
//
//Probes readProbes(const json& j)
//{
//	// TODO
//}
//
mfem::Mesh readMesh(const json& j)
{
	// TODO
	return mfem::Mesh("");
}

Model readModel(const json& j)
{
	auto modelJSON{ j.find("model") };
	if (modelJSON == j.end()) {
		throw std::runtime_error("Can not find \"model\" label.");
	}

	auto materialsJSON{ modelJSON->find("materials") };
	if (materialsJSON == modelJSON->end()) {
		throw std::runtime_error("Can not find \"materials\" label.");
	}
	
	auto pm{ semba::parsers::JSON::readMaterials(*materialsJSON) };
	
	AttributeToMaterial attToMat;
	for (const auto& mat : pm) {
		// TODO
	}

	AttributeToBoundary attToBdr;
	for (const auto& mat : pm) {
		// TODO
	}

	AttributeToBoundary attToInteriorConditions;
	for (const auto& mat : pm) {
		// TODO
	}


	return { readMesh(*modelJSON), attToMat, attToBdr, attToInteriorConditions };
}

Problem OpensembaAdapter::readProblem() const
{
	Problem r;
	auto model{ readModel(json_) };
	//auto sources{ readSources(json_)};
	//auto probes{ readProbes(json_) };
	

	return r;
}

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