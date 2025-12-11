#include "driver.h"
#include "string"

#include <numeric>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>

namespace maxwell::driver {

struct DSU { //Disjoint Set Union
    std::vector<int> p, r;
    explicit DSU(int n): p(n), r(n,0) { std::iota(p.begin(), p.end(), 0); }
    int find(int x){ return p[x]==x ? x : p[x]=find(p[x]); }
    void unite(int a, int b){
        a = find(a); b = find(b);
        if (a == b) return;
        if (r[a] < r[b]) std::swap(a,b);
        p[b] = a;
        if (r[a] == r[b]) r[a]++;
    }
};

std::vector<std::pair<int,int>> buildTwoElementPairsByTagToSort(mfem::Mesh& mesh, mfem::Array<int> tags)
{
	std::vector<std::pair<int,int>> res;
	for (auto t = 0; t < tags.Size(); t++){
		for(auto b = 0; b < mesh.GetNBE(); b++){
			if (mesh.GetBdrAttribute(b) == tags[t]){
				auto f_trans = mesh.GetInternalBdrFaceTransformations(b);
				if (auto f_trans = mesh.GetInternalBdrFaceTransformations(b)) {
    				res.emplace_back(f_trans->Elem1No, f_trans->Elem2No);
				}
			}
		}
	}
	return res;
}

static std::vector<std::pair<int,int>>
gatherConstraintPairs(mfem::Mesh& mesh,
                      const mfem::Array<int>& tfsf_tags,
                      const mfem::Array<int>& sbc_tags)
{
    std::vector<std::pair<int,int>> pairs;
    pairs.reserve(64);

    if (tfsf_tags.Size() > 0) {
        auto tfsf_pairs = buildTwoElementPairsByTagToSort(mesh, tfsf_tags);
        pairs.insert(pairs.end(), tfsf_pairs.begin(), tfsf_pairs.end());
    }
    if (sbc_tags.Size() > 0) {
        auto sbc_pairs = buildTwoElementPairsByTagToSort(mesh, sbc_tags);
        pairs.insert(pairs.end(), sbc_pairs.begin(), sbc_pairs.end());
    }
    return pairs;
}

DSU buildDSU(int NE, const std::vector<std::pair<int,int>>& pairs)
{
    DSU dsu(NE);
    for (const auto& pr : pairs) {
        const int a = pr.first, b = pr.second;
        if (0 <= a && a < NE && 0 <= b && b < NE) dsu.unite(a, b);
    }
    return dsu;
}

std::unordered_map<int, std::vector<int>> componentsFromDSU(int NE, DSU& dsu)
{
    std::unordered_map<int, std::vector<int>> comp;
    comp.reserve(NE);
    for (int e = 0; e < NE; ++e) comp[dsu.find(e)].push_back(e);
    return comp;
}

std::vector<long long> computeLoad(int P, int NE, const int* partitioning)
{
    std::vector<long long> load(P, 0);
    for (int e = 0; e < NE; ++e) {
        const int r = partitioning[e];
        if (0 <= r && r < P) load[r]++;
    }
    return load;
}

int chooseRankForComponent(const std::vector<int>& elems, int P, const int* partitioning, const std::vector<long long>& load)
{
    std::vector<int> hist(P, 0);
    for (int e : elems) {
        const int r = partitioning[e];
        if (0 <= r && r < P) hist[r]++;
    }

    int best_rank = -1, best_count = -1;
    for (int r = 0; r < P; ++r) {
        if (hist[r] > best_count) { best_count = hist[r]; best_rank = r; }
    }

    if (best_rank < 0) {
        best_rank = 0;
        for (int r = 1; r < P; ++r)
            if (load[r] < load[best_rank]) best_rank = r;
        return best_rank;
    }

    long long best_load = std::numeric_limits<long long>::max();
    int tie_pick = best_rank;
    for (int r = 0; r < P; ++r) {
        if (hist[r] == best_count && load[r] < best_load) {
            best_load = load[r];
            tie_pick = r;
        }
    }
    return tie_pick;
}

void assignComponent(const std::vector<int>& elems,
                int rank,
                int* partitioning,
                std::vector<long long>& load)
{
    for (int e : elems) partitioning[e] = rank;
    load[rank] += static_cast<long long>(elems.size());
}

void applyPairwiseConstraintsPartitioning(mfem::Mesh& mesh,
                                          int* partitioning,
                                          const mfem::Array<int>& tfsf_tags,
                                          const mfem::Array<int>& sbc_tags)
{
    const int NE = mesh.GetNE();
    const int P  = Mpi::WorldSize();
    if (NE == 0 || P <= 0) return;

    auto pairs = gatherConstraintPairs(mesh, tfsf_tags, sbc_tags);
    if (pairs.empty()) return;

    auto dsu = buildDSU(NE, pairs);
    auto comp = componentsFromDSU(NE, dsu);
    auto load = computeLoad(P, NE, partitioning);

    for (auto& kv : comp) {
        const auto& elems = kv.second;
        if (elems.size() <= 1) continue;
        const int rank = chooseRankForComponent(elems, P, partitioning, load);
        assignComponent(elems, rank, partitioning, load);
    }
}

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

const Direction assignFieldPol(const std::string& direction)
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

std::unique_ptr<InitialField> buildSphericalBesselJ6InitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0, 0.0, 1.0 }))
{
	Sources res;
	Source::Position center = Source::Position({ 0.0, 0.0, 0.0 });
	return std::make_unique<InitialField>(SphericalBesselJ6(), ft, p, center);
}

std::unique_ptr<TotalField> buildGaussianPlanewave(
	double spread,
	const Source::Position mean,
	const Source::Polarization& pol,
	const Source::Propagation& dir,
	const FieldType ft = FieldType::E
)
{
	Position projMean(3);
	projMean = 0.0;
	for (auto v = 0; v < mean.Size(); v++) {
		projMean[v] = mean[v];
	}
	Gaussian gauss{ spread, mfem::Vector({projMean * dir / dir.Norml2()})};
	Planewave pw(gauss, pol, dir, ft);
	return std::make_unique<TotalField>(pw);
}

std::unique_ptr<TotalField> buildDerivGaussDipole(
	const double length, 
	const double gaussianSpread, 
	const double gaussMean
) 
{
	DerivGaussDipole dip(length, gaussianSpread, gaussMean);
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
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6_2D") {
				res.add(buildBesselJ6InitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]))
				);
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6_3D") {
				res.add(buildSphericalBesselJ6InitialField(
					assignFieldType(case_data["sources"][s]["field_type"]),
					assemble3DVector(case_data["sources"][s]["polarization"]))
				);
			}
		}
		else if (case_data["sources"][s]["type"] == "planewave") {
			res.add(buildGaussianPlanewave(
				case_data["sources"][s]["magnitude"]["spread"],
				assemble3DVector(case_data["sources"][s]["magnitude"]["mean"]),
				assemble3DVector(case_data["sources"][s]["polarization"]),
				assemble3DVector(case_data["sources"][s]["propagation"]),
				FieldType::E)
			);
		}
		else if (case_data["sources"][s]["type"] == "dipole") {
			res.add(buildDerivGaussDipole(
				case_data["sources"][s]["magnitude"]["length"],
				case_data["sources"][s]["magnitude"]["spread"],
				case_data["sources"][s]["magnitude"]["mean"])
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

		if (case_data["solver_options"].contains("upwind_alpha")) {
			res.setUpwindAlpha(case_data["solver_options"]["upwind_alpha"]);
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
		
		if (case_data["solver_options"].contains("export_operator")) {
			res.setExportEO(case_data["solver_options"]["export_operator"]);
		}

		if (case_data["solver_options"].contains("evolution_operator")) {
			if (case_data["solver_options"]["evolution_operator"] == "maxwell") {
				res.setEvolutionOperator(EvolutionOperatorType::Maxwell);
			}
			else if (case_data["solver_options"]["evolution_operator"] == "global") {
				res.setEvolutionOperator(EvolutionOperatorType::Global);
			}
			else if (case_data["solver_options"]["evolution_operator"] == "hesthaven") {
				res.setEvolutionOperator(EvolutionOperatorType::Hesthaven);
			}
			else {
				throw std::runtime_error("Wrong type of Evolution Operator defined, please choose: 'maxwell', 'global' or 'hesthaven'. If none explicitly stated, it will default to 'global'.");
			}
		}

		if (case_data["solver_options"].contains("basis_type")) {
			res.setBasisType(case_data["solver_options"]["basis_type"]);
		}

		if (case_data["solver_options"].contains("tfsf_final_time")){
			res.setTFSFCutoffTime(case_data["solver_options"]["tfsf_final_time"]);
		}

		if (case_data["solver_options"].contains("ode_type")){
			res.setODEType(case_data["solver_options"]["ode_type"]);
		}

	}

	for (int b = 0; b < case_data["model"]["boundaries"].size(); ++b) {
        if (case_data["model"]["boundaries"][b].contains("type") && 
			case_data["model"]["boundaries"][b]["type"] == "SGBC"){ 
			if (!case_data["model"]["boundaries"][b].contains("material")){
				throw std::runtime_error("SGBC Material defined without material properties. Verify .json parameters.");
			} else {
				for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
					double rel_eps, rel_mu, sigma;
					if (!case_data["model"]["boundaries"][b]["material"].contains("relative_permittivity")){
						std::cout << "SGBC Material defined without 'relative_permittivity' parameter, assuming vacuum." << std::endl;
						rel_eps = 1.0;
					} 
					else{
						rel_eps = case_data["model"]["boundaries"][b]["material"]["relative_permittivity"];
					}
					if (!case_data["model"]["boundaries"][b]["material"].contains("relative_permeability")){
						std::cout << "SGBC Material defined without 'relative_permeability' parameter, assuming vacuum." << std::endl;
						rel_mu = 1.0;
					} 
					else{
						rel_mu = case_data["model"]["boundaries"][b]["material"]["relative_permeability"];
					}
					if (!case_data["model"]["boundaries"][b]["material"].contains("bulk_conductivity")){
						throw std::runtime_error("SGBC Material defined without 'bulk_conductivity parameter. Verify .json parameters.");
					} 
					else {
						sigma = case_data["model"]["boundaries"][b]["material"]["bulk_conductivity"]; // sigma_solver = sigma_si * Z0;
					}
					Material mat(rel_eps, rel_mu, sigma);
					SGBCProperties props(mat);
					for (auto t = 0; t < case_data["model"]["boundaries"][b]["tags"].size(); t++){
						props.geom_tags.emplace_back(case_data["model"]["boundaries"][b]["tags"][t]);
					}
					if (case_data["model"]["boundaries"][b]["material"].contains("num_of_segments")){
						props.num_of_segments = int(case_data["model"]["boundaries"][b]["material"]["num_of_segments"]);
					}
					if (case_data["model"]["boundaries"][b]["material"].contains("order")){
						props.order = int(case_data["model"]["boundaries"][b]["material"]["order"]);
					}
					if (case_data["model"]["boundaries"][b]["material"].contains("material_width")){
						props.material_width = double(case_data["model"]["boundaries"][b]["material"]["material_width"]);
					} 
					res.sgbc_props.emplace_back(props);
				}
			}
		}
	}

	return res;
}


Probes buildProbes(const json& case_data)
{

	Probes probes;
	size_t pp_count = 0;
	size_t fp_count = 0;

	if (case_data.contains("probes")){
		if (case_data["probes"].contains("exporter")) {
			ExporterProbe exporter_probe;
			if (case_data["probes"]["exporter"].contains("name")) {
				exporter_probe.name = case_data["probes"]["exporter"]["name"];
			}
			else{
				exporter_probe.name = case_data["model"]["filename"];
			}
			if (case_data["probes"]["exporter"].contains("steps")) {
				exporter_probe.visSteps = case_data["probes"]["exporter"]["steps"];
			}
			probes.exporterProbes.push_back(exporter_probe);
		}

		if (case_data["probes"].contains("point")) {
			for (int p = 0; p < case_data["probes"]["point"].size(); p++) {
				bool write = false;
				if (case_data["probes"]["point"][p].contains("write")){
					write = case_data["probes"]["point"][p]["write"];
				}
				PointProbe point_probe(
					assembleVector(case_data["probes"]["point"][p]["position"]),
					write
				);
				point_probe.setProbeID(pp_count);
				pp_count++;
				probes.pointProbes.push_back(point_probe);
			}
		}

			if (case_data["probes"].contains("field")) {
			for (int p = 0; p < case_data["probes"]["field"].size(); p++) {
				bool write = false;
				if (case_data["probes"]["field"][p].contains("write")){
					write = true;
				}
				FieldProbe field_probe(
					assignFieldType(case_data["probes"]["field"][p]["field_type"]),
					assignFieldPol(case_data["probes"]["field"][p]["polarization"]),
					assembleVector(case_data["probes"]["field"][p]["position"]),
					write
				);
				field_probe.setProbeID(fp_count);
				fp_count++;
				probes.fieldProbes.push_back(field_probe);
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
				if (case_data["probes"]["farfield"][p].contains("steps")) {
					probe.expSteps = case_data["probes"]["farfield"][p]["steps"];
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

		if (case_data["probes"].contains("domain_snapshot")){
			DomainSnapshotProbe probe;
			if (case_data["probes"]["domain_snapshot"].contains("name")) {
				probe.name = case_data["probes"]["domain_snapshot"]["name"];
			}
			else{
				probe.name = case_data["model"]["filename"];
			}
			if (case_data["probes"]["domain_snapshot"].contains("steps")) {
				probe.expSteps = case_data["probes"]["domain_snapshot"]["steps"];
			}
			probes.domainSnapshotProbes.push_back(probe);
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
	else if (bdr_cond == "SGBC") {
		return BdrCond::SGBC;
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
				if (case_data["model"]["boundaries"][b]["material"].contains("relative_permittivity")) {
					eps = case_data["model"]["boundaries"][b]["material"]["relative_permittivity"];
				}
				if (case_data["model"]["boundaries"][b]["material"].contains("relative_permeability")) {
					mu = case_data["model"]["boundaries"][b]["material"]["relative_permeability"];
				}
				if (case_data["model"]["boundaries"][b]["material"].contains("bulk_conductivity")) {
					sigma = case_data["model"]["boundaries"][b]["material"]["bulk_conductivity"];
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
		std::map<GeomTag, isInterior> geomTag2IsInterior;
	};

	checkBoundaryInputProperties(case_data);
	auto face2BdrEl{ mesh.GetFaceToBdrElMap() };

	geomTag2Info gt2i;

	for (auto b = 0; b < case_data["model"]["boundaries"].size(); b++) {
		for (auto a = 0; a < case_data["model"]["boundaries"][b]["tags"].size(); a++) {
			gt2i.geomTag2BdrCond.emplace(case_data["model"]["boundaries"][b]["tags"][a], assignBdrCond(case_data["model"]["boundaries"][b]["type"]));
			gt2i.geomTag2IsInterior.emplace(case_data["model"]["boundaries"][b]["tags"][a], false);
			for (auto f = 0; f < mesh.GetNumFaces(); f++) {
				if (face2BdrEl[f] != -1) {
					if (mesh.GetBdrAttribute(face2BdrEl[f]) == case_data["model"]["boundaries"][b]["tags"][a] 
						&& mesh.FaceIsInterior(f)) {
						gt2i.geomTag2IsInterior[case_data["model"]["boundaries"][b]["tags"][a]] = true;
						break;
					}
				}
			}
		}
	}

	GeomTagToBoundaryInfo res;
	for (auto [att, isInt] : gt2i.geomTag2IsInterior) {
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

Array<int> getTFSFTags(const json& case_data)
{
    Array<int> res;
    if (!case_data.contains("sources")) return res;

    for (int s = 0; s < case_data["sources"].size(); ++s) {
        if (case_data["sources"][s].contains("type") && case_data["sources"][s]["type"] == "planewave") {
            for (int t = 0; t < case_data["sources"][s]["tags"].size(); ++t) {
                res.Append(case_data["sources"][s]["tags"][t]);
            }
        }
    }
    return res;
}

Array<int> getSGBCTags(const json& case_data)
{
    Array<int> res;
    if (!case_data.contains("model") ||
        !case_data["model"].contains("boundaries")) {
        return res;
    }

    for (int b = 0; b < case_data["model"]["boundaries"].size(); ++b) {
        if (case_data["model"]["boundaries"][b].contains("type") && case_data["model"]["boundaries"][b]["type"] == "SGBC") {
            for (int t = 0; t < case_data["model"]["boundaries"][b]["tags"].size(); ++t) {
                res.Append(case_data["model"]["boundaries"][b]["tags"][t]);
            }
        }
    }
    return res;
}

Model buildModel(const json& case_data, const std::string& case_path, const bool isTest)
{
    mfem::Mesh mesh;
    if (isTest) {
        mesh = assembleMesh(assembleMeshString(case_data["model"]["filename"]));
    } else {
        mesh = assembleMesh(assembleLauncherMeshString(case_data["model"]["filename"], case_path));
    }

    auto att_to_material{ assembleAttributeToMaterial(case_data, mesh) };
    auto att_to_bdr_info{ assembleAttributeToBoundary(case_data, mesh) };

    if (!att_to_bdr_info.gt2b.empty() && !att_to_material.gt2m.empty()) {
        if (isTest) {
            mesh = assembleMesh(assembleMeshString(case_data["model"]["filename"]));
        } else {
            mesh = assembleMesh(assembleLauncherMeshString(case_data["model"]["filename"], case_path));
        }
    }

    if (case_data["model"].contains("refinement")) {
        auto ref_levels = int(case_data["model"]["refinement"]);
        for (auto r = 0; r < ref_levels; r++) {
            mesh.UniformRefinement();
        }
    }

    int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
    mfem::Array<int> tfsf_tags = getTFSFTags(case_data);
    mfem::Array<int> sbc_tags  = getSGBCTags(case_data);

	if (tfsf_tags.Size() || sbc_tags.Size()){
    	applyPairwiseConstraintsPartitioning(mesh, partitioning, tfsf_tags, sbc_tags);
	}

    Model res(mesh, att_to_material, att_to_bdr_info, partitioning);
    std::string filename = case_data["model"]["filename"];

    auto ends_with = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() &&
               str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (ends_with(filename, ".msh")) {
        res.meshName_ = filename.substr(0, filename.size() - 4);
    } else if (ends_with(filename, ".mesh")) {
        res.meshName_ = filename.substr(0, filename.size() - 5);
    } else {
        throw std::runtime_error("File format for mesh must be name.msh or name.mesh");
    }

    return res;
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

void postProcessInformation(const json& case_data, maxwell::Model& model, maxwell::SolverOptions& solverOpts) 
{
	for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		mfem::Array<int> tfsf_tags;
		if (case_data["sources"][s]["type"] == "planewave" || case_data["sources"][s]["type"] == "dipole") {
			for (auto t{ 0 }; t < case_data["sources"][s]["tags"].size(); t++) {
				tfsf_tags.Append(case_data["sources"][s]["tags"][t]);
			}
			auto tfsf_atts_present_in_partition_marker{ model.getMarker(maxwell::BdrCond::TotalFieldIn, true) };
			tfsf_atts_present_in_partition_marker.SetSize(model.getConstMesh().bdr_attributes.Max());
			tfsf_atts_present_in_partition_marker = 0;
			for (auto t = 0; t < tfsf_tags.Size(); t++){
				for (auto b = 0; b < model.getConstMesh().GetNBE(); b++){	
					if (model.getMesh().GetBdrAttribute(b) == tfsf_tags[t]){
						tfsf_atts_present_in_partition_marker[model.getMesh().GetBdrAttribute(b) - 1] = 1;
					}
				}
			}
			if (tfsf_atts_present_in_partition_marker.Sum() != 0){
				model.getTotalFieldScatteredFieldToMarker().insert(std::make_pair(maxwell::BdrCond::TotalFieldIn, tfsf_atts_present_in_partition_marker));
			}
		}
	}

	if (model.getBoundaryToMarker().find(BdrCond::SMA) != model.getBoundaryToMarker().end() && solverOpts.evolution.alpha == 0.0 && solverOpts.evolution.op == EvolutionOperatorType::Hesthaven) {
		throw std::runtime_error("Centered SMA with Hesthaven Evolution Operator not supported yet.");
	}

	model.getMesh().SetAttributes();
}

maxwell::Solver buildSolver(const json& case_data, const std::string& case_path, const bool isTest)
{
	
	maxwell::SolverOptions solverOpts{ buildSolverOptions(case_data) };
	maxwell::Sources sources{ buildSources(case_data) };
	maxwell::Probes probes{ buildProbes(case_data) };
	maxwell::Model model{ buildModel(case_data, case_path, isTest) };

	postProcessInformation(case_data, model, solverOpts);

	return maxwell::Solver(model, probes, sources, solverOpts);
}

}
