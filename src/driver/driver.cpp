#include "driver.h"

namespace maxwell::driver {

inline void checkIfThrows(bool condition, const std::string& msg)
{
	if (!condition) {
		throw std::runtime_error(msg.c_str());
	}
}

const FieldType assignFieldType(const std::string& field_type)
{
	if (field_type == "electric") {
		return FieldType::E;
	}
	else if (field_type == "magnetic") {
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

std::unique_ptr<InitialField> buildGaussianInitialField(
	const FieldType& ft = E,
	const double spread = 0.1,
	const mfem::Vector& center_ = mfem::Vector({ 0.5 }),
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const int dimension = 1)
{
	mfem::Vector gaussianCenter(dimension);
	gaussianCenter = 0.0;

	Gaussian gauss(spread, gaussianCenter, dimension);
	return std::make_unique<InitialField>(gauss, ft, p, center_);
}

std::unique_ptr<InitialField> buildResonantModeInitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const std::vector<std::size_t>& modes = { 1 })
{
	Sources res;
	Source::Position center((int)modes.size());
	center = 0.0;
	return std::make_unique<InitialField>(SinusoidalMode{ modes }, ft, p, center);
}

std::unique_ptr<InitialField> buildBesselJ6InitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0, 0.0, 1.0 }))
{
	Sources res;
	Source::Position center = Source::Position({ 0.0, 0.0, 0.0 });
	return std::make_unique<InitialField>(BesselJ6(), ft, p, center);
}

std::unique_ptr<TotalField> buildGaussianPlanewave(
	double spread,
	double delay,
	const Source::Polarization& pol,
	const Source::Propagation& dir,
	const FieldType ft = FieldType::E
)
{
	Gaussian gauss{ spread, mfem::Vector({-delay}) };
	Planewave pw(gauss, pol, dir, ft);
	return std::make_unique<TotalField>(pw);
}

std::unique_ptr<TotalField> buildDerivGaussDipole(
	const double length, 
	const double gaussianSpread, 
	const double gaussDelay
) 
{
	DerivGaussDipole dip(length, gaussianSpread, gaussDelay);
	return std::make_unique<TotalField>(dip);
}

Sources buildSources(const json& case_data)
{
	Sources res;
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] == "initial") {
			if (case_data["sources"][s]["magnitude"]["type"] == "gaussian") {
				res.add(buildGaussianInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					case_data["sources"][s]["magnitude"]["spread"],
					assembleCenterVector(case_data["sources"][s]["center"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["dimension"])
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "resonant") {
				res.add(buildResonantModeInitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]),
					case_data["sources"][s]["magnitude"]["modes"])
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6") {
				res.add(buildBesselJ6InitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]))
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "planewave") {
			res.add(buildGaussianPlanewave(
				case_data["sources"][s]["magnitude"]["spread"],
				case_data["sources"][s]["magnitude"]["delay"],
				assemble3DVector(case_data["sources"][s]["polarization"]),
				assemble3DVector(case_data["sources"][s]["propagation"]),
				FieldType::E)
			);
		}
		else if (case_data["sources"][s]["type"] == "dipole") {
			res.add(buildDerivGaussDipole(
				case_data["sources"][s]["magnitude"]["length"],
				case_data["sources"][s]["magnitude"]["spread"],
				case_data["sources"][s]["magnitude"]["delay"])
			);
		}
		else {
			throw std::runtime_error("Unknown source type in Json.");
		}
	}
	return res;
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

		if (case_data["solver_options"].contains("high_order_mesh")) {
			res.highOrderMesh = true;
		}

		if (case_data["solver_options"].contains("hesthaven_operator")) {
			res.setHesthavenOperator(case_data["solver_options"]["hesthaven_operator"]);
		}

		if (case_data["solver_options"].contains("global_operator")) {
			res.setGlobalOperator(case_data["solver_options"]["global_operator"]);
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
		if (case_data["probes"]["exporter"].contains("expSteps")) {
			exporter_probe.visSteps = case_data["probes"]["exporter"]["expSteps"];
		}
		probes.exporterProbes.push_back(exporter_probe);
	}

	if (case_data["probes"].contains("point")) {
		for (int p = 0; p < case_data["probes"]["point"].size(); p++) {
			PointProbe field_probe(
				assembleVector(case_data["probes"]["point"][p]["position"])
			);
			probes.pointProbes.push_back(field_probe);
		}
	}

	if (case_data["probes"].contains("farfield")) {
		for (int p{ 0 }; p < case_data["probes"]["farfield"].size(); p++) {
			NearFieldProbe probe;
			if (case_data["probes"]["farfield"][p].contains("name")) {
				probe.name = case_data["probes"]["farfield"][p]["name"];
			}
			if (case_data["probes"]["farfield"][p].contains("export_path")) {
				probe.exportPath = case_data["probes"]["farfield"][p]["export_path"];
			}
			if (case_data["probes"]["farfield"][p].contains("export_steps")) {
				probe.expSteps = case_data["probes"]["farfield"][p]["export_steps"];
			}
			if (case_data["probes"]["farfield"][p].contains("tags")) {
				std::vector<int> tags;
				for (int t{ 0 }; t < case_data["probes"]["farfield"][p]["tags"].size(); t++) {
					tags.push_back(case_data["probes"]["farfield"][p]["tags"][t]);
				}
				probe.tags = tags;
			}
			else {
				throw std::runtime_error("Tags have not been defined in farfield probe.");
			}
			probes.nearFieldProbes.push_back(probe);
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

std::string assembleLauncherMeshString(const std::string& mesh_name, const std::string& case_path)
{
	std::string path{ case_path };
	std::string s_json = ".json";
	std::string s_msh  = ".msh";
	std::string s_mesh = ".mesh";

	std::string::size_type input_json = path.find(s_json);
	std::string::size_type input_msh = mesh_name.find(s_msh);
	std::string::size_type input_mesh = mesh_name.find(s_mesh);

	path.erase(input_json, s_json.length());
	if (input_msh != std::string::npos)
		path.append(s_msh);
	if (input_mesh != std::string::npos)
		path.append(s_mesh);

	return path;
}

void checkIfAttributesArePresent(const Mesh& mesh, const GeomTagToMaterialInfo& info)
{

	for (auto [att, v] : info.gt2m) {
		checkIfThrows(
			mesh.attributes.Find(att) != -1,
			std::string("There is no attribute") + std::to_string(att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}

	for (auto [bdr_att, v] : info.gt2bm) {
		checkIfThrows(
			mesh.bdr_attributes.Find(bdr_att) != -1,
			std::string("There is no bdr_attribute") + std::to_string(bdr_att) +
			" defined in the mesh, but it is defined in the JSON."
		);
	}
}

GeomTagToMaterialInfo assembleAttributeToMaterial(const json& case_data, const mfem::Mesh& mesh)
{
	GeomTagToMaterialInfo res{};

	checkIfThrows(case_data.contains("model"), "JSON data does not include 'model'.");
	checkIfThrows(case_data["model"].contains("materials"), "JSON data does not include 'materials'.");

	for (auto m = 0; m < case_data["model"]["materials"].size(); m++) {
		for (auto t = 0; t < case_data["model"]["materials"][m]["tags"].size(); t++) {
			double eps{ 1.0 }, mu{ 1.0 }, sigma{ 0.0 };
			if (case_data["model"]["materials"][m].contains("relative_permittivity")) {
				eps = case_data["model"]["materials"][m]["relative_permittivity"];
			}
			if (case_data["model"]["materials"][m].contains("relative_permeability")) {
				mu = case_data["model"]["materials"][m]["relative_permeability"];
			}
			if (case_data["model"]["materials"][m].contains("bulk_conductivity")) {
				sigma = case_data["model"]["materials"][m]["bulk_conductivity"];
			}
			res.gt2m.emplace(std::make_pair(case_data["model"]["materials"][m]["tags"][t], Material(eps, mu, sigma)));
		}
	}

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			if (case_data["model"]["boundaries"][b].contains("material")) {
				double eps{ 1.0 }, mu{ 1.0 }, sigma{ 0.0 };
				if (case_data["model"]["boundaries"][b].contains("relative_permittivity")) {
					eps = case_data["model"]["boundaries"][b]["relative_permittivity"];
				}
				if (case_data["model"]["boundaries"][b].contains("relative_permeability")) {
					mu = case_data["model"]["boundaries"][b]["relative_permeability"];
				}
				if (case_data["model"]["boundaries"][b].contains("bulk_conductivity")) {
					sigma = case_data["model"]["boundaries"][b]["bulk_conductivity"];
				}
				res.gt2bm.emplace(case_data["model"]["boundaries"][b]["tags"][a], Material(eps, mu, sigma));
			}
		}
	}

	checkIfAttributesArePresent(mesh, res);

	return res;
}

void checkBoundaryInputProperties(const json& case_data)
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
}

GeomTagToBoundaryInfo assembleAttributeToBoundary(const json& case_data, const mfem::Mesh& mesh)
{
	using isInterior = bool;

	struct geomTag2Info {
		std::map<GeomTag, BdrCond> geomTag2BdrCond;
		std::map<GeomTag, isInterior> geomTag2Interior;
	};

	checkBoundaryInputProperties(case_data);
	auto face2BdrEl{ mesh.GetFaceToBdrElMap() };

	geomTag2Info gt2i;

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			gt2i.geomTag2BdrCond.emplace(case_data["model"]["boundaries"][b]["tags"][a], assignBdrCond(case_data["model"]["boundaries"][b]["type"]));
			gt2i.geomTag2Interior.emplace(case_data["model"]["boundaries"][b]["tags"][a], false);
			for (auto f = 0; f < mesh.GetNumFaces(); f++) {
				if (face2BdrEl[f] != -1) {
					if (mesh.GetBdrAttribute(face2BdrEl[f]) == case_data["model"]["boundaries"][b]["tags"][a] 
						&& mesh.FaceIsInterior(f)) {
						gt2i.geomTag2Interior[case_data["model"]["boundaries"][b]["tags"][a]] = true;
						break;
					}
				}
			}
		}
	}

	GeomTagToBoundaryInfo res;
	for (auto [att, isInt] : gt2i.geomTag2Interior) {
		switch (isInt) {
		case false:
			res.gt2b.emplace(att, gt2i.geomTag2BdrCond[att]);
			break;
		case true:
			res.gt2ib.emplace(att, gt2i.geomTag2BdrCond[att]);
			break;
		}
	}

	return res;
}
mfem::Mesh assembleMesh(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFile(mesh_string, 1, 0, false);
}

mfem::Mesh assembleMeshNoFix(const std::string& mesh_string)
{
	return mfem::Mesh::LoadFromFileNoBdrFix(mesh_string, 1, 0, false);
}

Model buildModel(const json& case_data, const std::string& case_path, const bool isTest)
{
	mfem::Mesh mesh;
	if (isTest) {
		mesh = assembleMesh(assembleMeshString(case_data["model"]["filename"]));
	}
	else {
		mesh = assembleMesh(assembleLauncherMeshString(case_data["model"]["filename"], case_path));
	}

	auto att_to_material{ assembleAttributeToMaterial(case_data, mesh) };
	auto att_to_bdr_info{ assembleAttributeToBoundary(case_data, mesh) };

	if (!att_to_bdr_info.gt2b.empty() && !att_to_material.gt2m.empty()) {
		if (isTest) {
			mesh = assembleMesh(assembleMeshString(case_data["model"]["filename"]));
		}
		else {
			mesh = assembleMesh(assembleLauncherMeshString(case_data["model"]["filename"], case_path));
		}
	}

	return Model(mesh, att_to_material, att_to_bdr_info);
}

json parseJSONfile(const std::string& case_name)
{
	std::ifstream test_file(case_name);
	return json::parse(test_file);
}

maxwell::Solver buildSolverJson(const std::string& case_name, const bool isTest)
{
	auto case_data = parseJSONfile(case_name);

	return buildSolver(case_data, case_name, isTest);
}

void postProcessInformation(const json& case_data, maxwell::Model& model) 
{
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		mfem::Array<int> tfsf_tags;
		if (case_data["sources"][s]["type"] == "planewave" || case_data["sources"][s]["type"] == "dipole") {
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

maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest)
{
	
	maxwell::Model model{ buildModel(case_data, case_path, isTest) };
	maxwell::Probes probes{ buildProbes(case_data) };
	maxwell::Sources sources{ buildSources(case_data) };
	maxwell::SolverOptions solverOpts{ buildSolverOptions(case_data) };

	postProcessInformation(case_data, model);

	return maxwell::Solver(model, probes, sources, solverOpts);
}

}