#include "CasesFunctions.h"
#include <regex>
#include <nlohmann/json.hpp>

namespace maxwell{

    using namespace mfem;
    using json = nlohmann::json;

json parseJSONfile(const std::string& case_name)
{
	std::ifstream test_file(case_name);
	return json::parse(test_file);
}

Mesh loadMeshFromFile(const std::string& mesh_path)
{
    Mesh res;
    res.LoadFromFile(mesh_path, 0, 0, 0);
    res.Finalize();
    return res;
}

GridFunction loadGridFunctionFromFile(const std::string& file_path, Mesh& mesh)
{
    std::ifstream in(file_path);
	GridFunction res(&mesh, in);
	return res;
}

int extractRankNumber(const std::string& folder_name) {
    static std::regex rank_regex(R"(rank_(\d+))");
    std::smatch match;
    if (std::regex_match(folder_name, match, rank_regex)) {
        return std::stoi(match[1].str());
    }
    throw std::runtime_error("Invalid rank folder name: " + folder_name);
}

const std::vector<Position> buildDoFPositions(const GridFunction& gf)
{
    auto fes = gf.FESpace();
    auto fec{ dynamic_cast<const L2_FECollection*>(fes->FEColl()) };
    auto vdimfes = FiniteElementSpace(fes->GetMesh(), fec, 3);
    GridFunction nodes(&vdimfes);
    fes->GetMesh()->GetNodes(nodes);
    auto dirSize{ nodes.Size() / 3 };
    std::vector<Position> res;
    res.resize(dirSize);
    for (auto i{ 0 }; i < dirSize; i++) {
        res[i] = Position({ nodes[i], nodes[i + dirSize], nodes[i + dirSize * 2] });
    }
    return res;
}

void L2SimDataCalculator::loadNodepos(const std::string& data_path)
{
    std::filesystem::path probes_dir = std::filesystem::path(data_path) / "DomainSnapshopProbes";
    if (!std::filesystem::exists(probes_dir) || !std::filesystem::is_directory(probes_dir)) {
        throw std::runtime_error("No DomainSnapshopProbes were generated in this case. Rerun case with one.");
    }

    for (const auto& entry : std::filesystem::directory_iterator(probes_dir)) {
        if (!entry.is_directory()) continue;

        std::string folder_name = entry.path().filename().string();

        if (folder_name == "meshes") continue;

        if (folder_name.rfind("rank_", 0) == 0) {
            int rank_no = extractRankNumber(folder_name);

            std::filesystem::path ex_file = entry.path() / "cycle_000000" / "Ex.gf";
            if (!std::filesystem::exists(ex_file)) {
                std::cerr << "Warning: Ex.gf missing in " << ex_file.parent_path() << "\n";
                continue;
            }

            if (rank_no >= (int)meshes_.size()) {
                throw std::runtime_error("Mesh not available for rank " + std::to_string(rank_no));
            }

            nodepos_[rank_no] = buildDoFPositions(loadGridFunctionFromFile(ex_file.string(), meshes_[rank_no]));
        }
    }
}

void L2SimDataCalculator::loadMeshes(const std::string& data_path)
{
    std::filesystem::path meshPath = std::filesystem::path(data_path) / "/DomainSnapshopProbes/meshes";
    if (!std::filesystem::exists(meshPath) || !std::filesystem::is_directory(meshPath)) {
        throw std::runtime_error("No DomainSnapshopProbes were generated in this case. Rerun case with one.");
    }

    int rank_no = 0;
    for (const auto& entry : std::filesystem::directory_iterator(meshPath)) {
        if (entry.is_regular_file()) {
            const auto& fname = entry.path().filename().string();
            if (fname.rfind("mesh_rank", 0) == 0) {
                meshes_[rank_no] = loadMeshFromFile(entry.path().string());
                rank_no++;
            }
        }
    }
}

FunctionType getFunctionTypeFromJson(const json& case_data)
{
    for (auto s{ 0 }; s < case_data["sources"].size(); s++) {
		if (case_data["sources"][s]["type"] == "initial") {
			if (case_data["sources"][s]["magnitude"]["type"] == "gaussian") {
                return FunctionType::Gaussian; //Currently unsupported.
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "resonant") {
                return FunctionType::Resonant;
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6_2D") {
                return FunctionType::BesselJ62D; //Currently unsupported.
			}
			else if (case_data["sources"][s]["magnitude"]["type"] == "besselj6_3D") {
                return FunctionType::BesselJ63D; //Currently unsupported.
            }
		}
		else if (case_data["sources"][s]["type"] == "planewave") {
            return FunctionType::Planewave; //Currently unsupported.
		}
		else if (case_data["sources"][s]["type"] == "dipole") {
            return FunctionType::Dipole; //Currently unsupported.
		}
		else {
			throw std::runtime_error("Unknown source type in Json.");
		}
	}
    throw std::runtime_error("No source detected in getFunctionTypeFromJson.");
}

std::unique_ptr<TimeFunction> buildFunctionByType(const json& case_data)
{
    if(getFunctionTypeFromJson(case_data) == FunctionType::Resonant){
        std::vector<int> modes = case_data["sources"][0]["magnitude"]["modes"];
        std::vector<double> box_size(modes.size()); //Assuming box size equal to 1.0 x 1.0 (x1.0 if 3D)
        for (auto v = 0; v < box_size.size(); v++){
            box_size[v] = 1.0;
        }
        return std::make_unique<TimeResonantSinusoidalMode>(case_data["sources"][0]["magnitude"]["modes"], box_size);
    }
    else{
        throw std::runtime_error("Currently unsupported FunctionType in buildFunctionByType.");
    }
}

void L2SimDataCalculator::initFunction(const std::string& json_path)
{
    auto case_data = parseJSONfile(json_path);
    function_ = buildFunctionByType(case_data);
}

L2SimDataCalculator::L2SimDataCalculator(const std::string& data_path, const std::string& json_file)
{
    loadMeshes(data_path);
    loadNodepos(data_path);
    
    initFunction(json_file);
    


}



}