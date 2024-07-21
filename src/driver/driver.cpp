#include "MaxwellAdapter.hpp"

inline void checkIfThrows(bool condition, const std::string& msg)
{
	if (!condition) {
		throw std::runtime_error(msg.c_str());
	}
}

const FieldType assignFieldType(const std::string& field_type)
{
	if (field_type == "E") {
		return FieldType::E;
	}
	else if (field_type == "H") {
		return FieldType::H;
	}
	else {
		throw std::runtime_error("Wrong Field Type in Point Probe assignation.");
	}
}

const Direction assignFieldSpatial(const std::string& direction)
{
	if (direction == "X") {
		return X;
	}
	else if (direction == "Y") {
		return Y;
	}
	else if (direction == "Z") {
		return Z;
	}
	else {
		throw std::runtime_error("Wrong Field Polarization in Point Probe assignation.");
	}
}

std::vector<double> assembleVector(const json& input)
{
	std::vector<double> res(input.size());
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
}

mfem::Vector assembleCenterVector(const json& source_center)
{
	mfem::Vector res(int(source_center.size()));
	for (int i = 0; i < source_center.size(); i++) {
		res[i] = source_center[i];
	}
	return res;
}

mfem::Vector assemble3DVector(const json& input)
{
	mfem::Vector res(3);
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
}

FieldType getFieldType(const std::string& ft)
{
	if (ft == "electric") {
		return FieldType::E;
	}
	else if (ft == "magnetic") {
		return FieldType::H;
	}
	else {
		throw std::runtime_error("The fieldtype written in the json is neither 'electric' nor 'magnetic'");
	}
}

Sources buildSources(const json& case_data)
{
	Sources res;
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] == "initial") {
			if (case_data["sources"][s]["magnitude"]["type"] == "gaussian") {
				return buildGaussianInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					case_data["sources"][s]["magnitude"]["spread"],
					assembleCenterVector(case_data["sources"][s]["center"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["dimension"]
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "resonant") {
				return buildResonantModeInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["magnitude"]["modes"]
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "totalField") {
			return buildGaussianPlanewave(
				case_data["sources"][s]["magnitude"]["spread"],
				case_data["sources"][s]["magnitude"]["delay"],
				assemble3DVector(case_data["sources"][s]["polarization"]),
				assemble3DVector(case_data["sources"][s]["propagation"]),
				getFieldType(case_data["sources"][s]["fieldtype"])
			);
		}
		else {
			throw std::runtime_error("Unknown source type in Json.");
		}
	}
}


SolverOptions buildSolverOptions(const json& case_data)
{
	SolverOptions res{};

	if (case_data.contains("solver_options")) {

		if (case_data["solver_options"].contains("solver_type")) {
			if (case_data["solver_options"]["solver_type"] == "centered")
				res.setCentered();
			else if (case_data["solver_options"]["solver_type"] == "upwind") {}
		}

		if (case_data["solver_options"].contains("time_step")) {
			res.setTimeStep(case_data["solver_options"]["time_step"]);
		}

		if (case_data["solver_options"].contains("final_time")) {
			res.setFinalTime(case_data["solver_options"]["final_time"]);
		}

		if (case_data["solver_options"].contains("cfl")) {
			res.setCFL(case_data["solver_options"]["cfl"]);
		}

		if (case_data["solver_options"].contains("order")) {
			res.setOrder(case_data["solver_options"]["order"]);
		}

		if (case_data["solver_options"].contains("spectral")) {
			res.setSpectralEO(case_data["solver_options"]["spectral"]);
		}

	}
	return res;
}


Probes buildProbes(const json& case_data)
{

	Probes probes;

	if (case_data["probes"].contains("exporter")) {
		ExporterProbe exporter_probe;
		exporter_probe.name = case_data["model"]["filename"];
		if (case_data["probes"]["exporter"].contains("steps")) {
			exporter_probe.visSteps = case_data["probes"]["exporter"]["steps"];
		}
		probes.exporterProbes.push_back(exporter_probe);
	}

	if (case_data["probes"].contains("field")) {
		for (int p = 0; p < case_data["probes"]["field"].size(); p++) {
			FieldProbe field_probe(
				assembleVector(case_data["probes"]["field"][p]["position"])
			);
			probes.fieldProbes.push_back(field_probe);
		}
	}

	if (case_data["probes"].contains("neartofarfield")) {
		for (int p{ 0 }; p < case_data["probes"]["neartofarfield"].size(); p++) {
			NearToFarFieldProbe probe;
			if (case_data["probes"]["neartofarfield"][p].contains("name")) {
				probe.name = case_data["probes"]["neartofarfield"][p]["name"];
			}
			if (case_data["probes"]["neartofarfield"][p].contains("steps")) {
				probe.steps = case_data["probes"]["neartofarfield"][p]["steps"];
			}
			if (case_data["probes"]["neartofarfield"][p].contains("tags")) {
				std::vector<int> tags;
				for (int t{ 0 }; t < case_data["probes"]["neartofarfield"][p]["tags"].size(); t++) {
					tags.push_back(case_data["probes"]["neartofarfield"][p]["tags"][t]);
				}
				probe.tags = tags;
			}
			else {
				throw std::runtime_error("Tags have not been defined in neartofarfield probe.");
			}
			probes.nearToFarFieldProbes.push_back(probe);
		}
	}

	return probes;
}


std::string dataFolder() { return "./testData/"; }
std::string maxwellInputsFolder() { return dataFolder() + "maxwellInputs/"; }

BdrCond assignBdrCond(const std::string& bdr_cond)
{
	if (bdr_cond == "PEC") {
		return BdrCond::PEC;
	}
	else if (bdr_cond == "PMC") {
		return BdrCond::PMC;
	}
	else if (bdr_cond == "SMA") {
		return BdrCond::SMA;
	}
	else {
		throw std::runtime_error(("The defined Boundary Type " + bdr_cond + " is incorrect.").c_str());
	}
}

std::string assembleMeshString(const std::string& filename)
{
	std::string folder_name{ filename };
	std::string s_msh = ".msh";
	std::string s_mesh = ".mesh";

	std::string::size_type input_msh = folder_name.find(s_msh);
	std::string::size_type input_mesh = folder_name.find(s_mesh);

	if (input_msh != std::string::npos)
		folder_name.erase(input_msh, s_msh.length());
	if (input_mesh != std::string::npos)
		folder_name.erase(input_mesh, s_mesh.length());

	return maxwellInputsFolder() + folder_name + "/" + filename;
}

GeomTagToMaterial assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh)
{
	GeomTagToMaterial res{};

	checkIfThrows(case_data.contains("model"), "JSON data does not include 'model'.");
	checkIfThrows(case_data["model"].contains("materials"), "JSON data does not include 'materials'.");

	auto show{ case_data["model"]["materials"] };

	for (auto m = 0; m < case_data["model"]["materials"].size(); m++) {
		for (auto t = 0; t < case_data["model"]["materials"][m]["tags"].size(); t++) {
			double eps{ 1.0 };
			double mu{ 1.0 };
			if (case_data["model"]["materials"][m].contains("relative_permittivity")) {
				eps = case_data["model"]["materials"][m]["relative_permittivity"];
			}
			if (case_data["model"]["materials"][m].contains("relative_permeability")) {
				mu = case_data["model"]["materials"][m]["relative_permeability"];
			}
			res.emplace(std::make_pair(case_data["model"]["materials"][m]["tags"][t], Material(eps, mu)));
		}
	}

	for (auto [att, v] : res) {
		checkIfThrows(
			mesh.attributes.Find(att) != -1,
			std::string("There is no attribute") + std::to_string(att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}

	//for (auto att{ 1 }; att < mesh.attributes.Size() + 1; att++) {
	//	checkIfThrows(
	//		!(res.find(att) == res.end()),
	//		std::string("There is no attribute") + std::to_string(att) +
	//		" defined in the JSON, but it is defined in the mesh."
	//	);
	//}

	return res;
}

BoundaryPair assembleAttributeToBoundary(const json& case_data, const mfem::Mesh& mesh)
{
	checkIfThrows(case_data["model"].contains("boundaries"),
		"JSON data does not include 'boundaries' in 'model'.");

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {

		checkIfThrows(case_data["model"]["boundaries"][b].contains("tags"),
			"Boundary " + std::to_string(b) + " does not have defined 'tags'.");

		checkIfThrows(!case_data["model"]["boundaries"][b]["tags"].empty(),
			"Boundary " + std::to_string(b) + " 'tags' are empty.");

		checkIfThrows(case_data["model"]["boundaries"][b].contains("type"),
			"Boundary " + std::to_string(b) + " does not have a defined 'type'.");
	}

	auto face2BdrEl{ mesh.GetFaceToBdrElMap() };

	std::map<GeomTag, BdrCond> geomTag2bdrCond;
	std::map<GeomTag, isInterior> geomTag2interior;
	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			geomTag2bdrCond.emplace(case_data["model"]["boundaries"][b]["tags"][a], assignBdrCond(case_data["model"]["boundaries"][b]["type"]));
			geomTag2interior.emplace(case_data["model"]["boundaries"][b]["tags"][a], false);
			for (auto f = 0; f < mesh.GetNumFaces(); f++) {
				if (face2BdrEl[f] != -1) {
					if (mesh.GetBdrAttribute(face2BdrEl[f]) == a && mesh.FaceIsInterior(f)) {
						geomTag2interior[a] = true;
						break;
					}
				}
			}
		}
	}

	//for (auto [geomTag, v] : geomTag2interior) {
	//	checkIfThrows(
	//		mesh.bdr_attributes.Find(geomTag) != -1,
	//		std::string("There is no boundary geometrical tag ") + std::to_string(geomTag) +
	//		" defined in the mesh, but it is defined in the JSON."
	//	);
	//}

	//for (auto geomTag{ 1 }; geomTag < mesh.bdr_attributes.Size() + 1; geomTag++) {
	//	checkIfThrows(
	//		!(geomTag2interior.find(geomTag) == geomTag2interior.end()),
	//		std::string("There is no  boundary geometrical tag ") + std::to_string(geomTag) +
	//		" defined in the JSON, but it is defined in the mesh."
	//	);
	//}

	GeomTagToBoundary geomTag2bdr{};
	GeomTagToInteriorConditions geomTag2intCond{};
	for (auto [att, isInt] : geomTag2interior) {
		switch (isInt) {
		case false:
			geomTag2bdr.emplace(att, geomTag2bdrCond[att]);
			break;
		case true:
			geomTag2intCond.emplace(att, geomTag2bdrCond[att]);
			break;
		}
	}

	return std::make_pair(geomTag2bdr, geomTag2intCond);
}

mfem::Mesh assembleMesh(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFile(mesh_string, 1, 0, true);
}

mfem::Mesh assembleMeshNoFix(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFileNoBdrFix(mesh_string, 1, 0, true);
}

Model buildModel(const json& case_data)
{
	auto mesh{ assembleMesh(assembleMeshString(case_data["model"]["filename"])) };

	auto att_to_material{ assembleAttributeToMaterial(case_data, mesh) };
	auto att_to_bdr_info{ assembleAttributeToBoundary(case_data, mesh) };

	if (!att_to_bdr_info.second.empty()) {
		mesh = assembleMeshNoFix(assembleMeshString(case_data["model"]["filename"]));
	}

	return Model(mesh, att_to_material, att_to_bdr_info.first, att_to_bdr_info.second);
}

json parseJSONfile(const std::string& case_name)
{
	auto file_name{ maxwellCase(case_name) };
	std::ifstream test_file(file_name);
	return json::parse(test_file);
}

maxwell::Solver buildSolver(const std::string& case_name)
{
	auto case_data{ parseJSONfile(case_name) };

	return buildSolver(case_data);
}

void postProcessInformation(const json& case_data, maxwell::Model& model) 
{
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		mfem::Array<int> tfsf_tags;
		if (case_data["sources"][s]["type"] == "totalField") {
			for (auto t{ 0 }; t < case_data["sources"][s]["tags"].size(); t++) {
				tfsf_tags.Append(case_data["sources"][s]["tags"][t]);
			}
		auto marker{ model.getMarker(maxwell::BdrCond::TotalFieldIn, true) };
		marker.SetSize(model.getConstMesh().bdr_attributes.Max());
		marker = 0;
		for (auto t : tfsf_tags) {
			marker[t - 1] = 1;
		}
		model.getTotalFieldScatteredFieldToMarker().insert(std::make_pair(maxwell::BdrCond::TotalFieldIn, marker));
		}
	}
}

maxwell::Solver buildSolver(const json& case_data)
{
	maxwell::Model model{ maxwell::buildModel(case_data) };
	maxwell::Probes probes{ maxwell::buildProbes(case_data) };
	maxwell::Sources sources{ maxwell::buildSources(case_data) };
	maxwell::SolverOptions solverOpts{ maxwell::buildSolverOptions(case_data) };

	postProcessInformation(case_data, model);

	return maxwell::Solver(model, probes, sources, solverOpts);
}
