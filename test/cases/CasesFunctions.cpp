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


Vector assemble3DVector(const json& input)
{
	Vector res(3);
	for (int i = 0; i < input.size(); i++) {
		res[i] = input[i];
	}
	return res;
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

int getExcitedDirectionsFromInitial(const Vector& pol)
{
    int res = 1;
    for (auto v = 0; v < pol.Size(); v++){
        if (pol[v] != 0.0) {
            res++;
        }
    }
    return res;
}

ExcitationCoeffs::ExcitationCoeffs(const std::string& json_path)
{
    auto case_data = parseJSONfile(json_path);
    auto ft = assignFieldType(case_data["sources"][0]["field_type"]);
    auto pol = assemble3DVector(case_data["sources"][0]["polarization"]);
    pol /= pol.Norml2();

    initFieldCompFactor();
    loadInitialPolarizationValues(ft, pol);
    loadExcitedDirectionValues(ft, pol);
}

void ExcitationCoeffs::initFieldCompFactor()
{
    FieldCompFactor[E][X] = 0.0; 
    FieldCompFactor[E][Y] = 0.0; 
    FieldCompFactor[E][Z] = 0.0; 
    FieldCompFactor[H][X] = 0.0; 
    FieldCompFactor[H][Y] = 0.0; 
    FieldCompFactor[H][Z] = 0.0; 
}

void ExcitationCoeffs::loadInitialPolarizationValues(const FieldType& ft, const Vector& pol)
{
    switch(ft){
        case E:
            if(pol[X] != 0.0){ FieldCompFactor[E][X] = pol[X];} else {FieldCompFactor[E][X] = 0.0;}
            if(pol[Y] != 0.0){ FieldCompFactor[E][Y] = pol[Y];} else {FieldCompFactor[E][Y] = 0.0;}
            if(pol[Z] != 0.0){ FieldCompFactor[E][Z] = pol[Z];} else {FieldCompFactor[E][Z] = 0.0;}
            break;
        case H:
            if(pol[X] != 0.0){ FieldCompFactor[H][X] = pol[X];} else {FieldCompFactor[H][X] = 0.0;}
            if(pol[Y] != 0.0){ FieldCompFactor[H][Y] = pol[Y];} else {FieldCompFactor[H][Y] = 0.0;}
            if(pol[Z] != 0.0){ FieldCompFactor[H][Z] = pol[Z];} else {FieldCompFactor[H][Z] = 0.0;}
            break;
    }
}

void normaliseExcitedVector(double& xcomp, double& ycomp, double& zcomp)
{
    auto normfactor = std::sqrt(std::pow(xcomp, 2.0) + std::pow(ycomp, 2.0) + std::pow(zcomp, 2.0));
    xcomp /= normfactor; 
    ycomp /= normfactor; 
    zcomp /= normfactor; 
}

void ExcitationCoeffs::loadExcitedDirectionValues(const FieldType& ft, const Vector& pol)
{
    auto numExcitedDirections = getExcitedDirectionsFromInitial(pol);
    if (numExcitedDirections == 2){
        switch(ft){
            case E:
                if (pol[X] != 0.0)     { FieldCompFactor[H][Y] = 1.0/std::sqrt(2.0); FieldCompFactor[H][Z] = 1.0/std::sqrt(2.0);}
                else if (pol[Y] != 0.0){ FieldCompFactor[H][X] = 1.0/std::sqrt(2.0); FieldCompFactor[H][Z] = 1.0/std::sqrt(2.0);}
                else                   { FieldCompFactor[H][X] = 1.0/std::sqrt(2.0); FieldCompFactor[H][Y] = 1.0/std::sqrt(2.0);}
                break;
            case H:
                if (pol[X] != 0.0)     { FieldCompFactor[E][Y] = 1.0/std::sqrt(2.0); FieldCompFactor[E][Z] = 1.0/std::sqrt(2.0);}
                else if (pol[Y] != 0.0){ FieldCompFactor[E][X] = 1.0/std::sqrt(2.0); FieldCompFactor[E][Z] = 1.0/std::sqrt(2.0);}
                else                   { FieldCompFactor[E][X] = 1.0/std::sqrt(2.0); FieldCompFactor[E][Y] = 1.0/std::sqrt(2.0);}
                break;
        }
    }
    else if (numExcitedDirections == 3){
        double normfactor = 0.0;
        switch(ft){
            case E:
                if (pol[X] != 0.0){ FieldCompFactor[H][Y] += pol[X]; FieldCompFactor[H][Z] += pol[X];}
                if (pol[Y] != 0.0){ FieldCompFactor[H][X] += pol[Y]; FieldCompFactor[H][Z] += pol[Y];}
                if (pol[Z] != 0.0){ FieldCompFactor[H][X] += pol[Z]; FieldCompFactor[H][Y] += pol[Z];}
                normaliseExcitedVector(FieldCompFactor[H][X], FieldCompFactor[H][Y], FieldCompFactor[H][Z]);
                break;
            case H:
                if (pol[X] != 0.0){ FieldCompFactor[E][Y] += pol[X]; FieldCompFactor[E][Z] += pol[X];}
                if (pol[Y] != 0.0){ FieldCompFactor[E][X] += pol[Y]; FieldCompFactor[E][Z] += pol[Y];}
                if (pol[Z] != 0.0){ FieldCompFactor[E][X] += pol[Z]; FieldCompFactor[E][Y] += pol[Z];}
                normaliseExcitedVector(FieldCompFactor[E][X], FieldCompFactor[E][Y], FieldCompFactor[E][Z]);
                break;
        }
    }
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

GridFunction getGridFunction(const std::string& grid_path, Mesh& mesh)
{
    std::ifstream in(grid_path);
	return GridFunction(&mesh, in);
}

double getTime(const std::filesystem::path& time_path)
{
    double res = 0.0;
    std::ifstream in(time_path);
    if (in >> res) {
        std::cout << "Simulation time = " << time_path << '\n';
    } else {
        std::cerr << "Failed to read a double from " << time_path << '\n';
    }
    return res;
}

int getGlobalNdofs(const std::map<Rank, std::vector<Position>>& nodepos)
{
    int res = 0;
    for (auto r = 0; r < nodepos.size(); r++)
    {
        res += nodepos.at(r).size();
    }
    return res;
}

L2SimDataCalculator::L2SimDataCalculator(const std::string& data_path, const std::string& json_path)
{
    loadMeshes(data_path);
    loadNodepos(data_path);
    initFunction(json_path);
    
    ExcitationCoeffs excCoeff(json_path);
    
    double l2diff = 0.0;
    int ndofs = getGlobalNdofs(nodepos_);

    for (auto r = 0; r < meshes_.size(); r++)
    {
        std::string rank_string("rank_" + std::to_string(r));
        std::filesystem::path rankPath = std::filesystem::path(data_path) / "/DomainSnapshopProbes/" / rank_string;

        for (auto& cycleEntry : std::filesystem::directory_iterator(rankPath)) {

            std::filesystem::path cyclePath = cycleEntry.path();

            std::filesystem::path ExPath = cyclePath / "Ex.gf";
            std::filesystem::path EyPath = cyclePath / "Ey.gf";
            std::filesystem::path EzPath = cyclePath / "Ez.gf";
            std::filesystem::path HxPath = cyclePath / "Hx.gf";
            std::filesystem::path HyPath = cyclePath / "Hy.gf";
            std::filesystem::path HzPath = cyclePath / "Hz.gf";
            std::filesystem::path timePath = cyclePath / "time.txt";

            auto Ex = getGridFunction(ExPath.string(), meshes_[r]);
            auto Ey = getGridFunction(EyPath.string(), meshes_[r]);
            auto Ez = getGridFunction(EzPath.string(), meshes_[r]);
            auto Hx = getGridFunction(HxPath.string(), meshes_[r]);
            auto Hy = getGridFunction(HyPath.string(), meshes_[r]);
            auto Hz = getGridFunction(HzPath.string(), meshes_[r]);
            double time = getTime(timePath);

            Vector analytic(nodepos_[r].size()); 

            for (auto v = 0; v < Ex.Size(); v++){
                analytic[v] = function_->eval(nodepos_[r][v], time);
            }

            l2diff += std::sqrt(Ex.Add(-1.0 * excCoeff.FieldCompFactor[E][X], analytic).Sum())
            + std::sqrt(Ey.Add(-1.0 * excCoeff.FieldCompFactor[E][Y], analytic).Sum())
            + std::sqrt(Ez.Add(-1.0 * excCoeff.FieldCompFactor[E][Z], analytic).Sum())
            + std::sqrt(Hx.Add(-1.0 * excCoeff.FieldCompFactor[H][X], analytic).Sum())
            + std::sqrt(Hy.Add(-1.0 * excCoeff.FieldCompFactor[H][Y], analytic).Sum())
            + std::sqrt(Hz.Add(-1.0 * excCoeff.FieldCompFactor[H][Z], analytic).Sum());

            l2diff /= double(ndofs);
        }
    }
}


}