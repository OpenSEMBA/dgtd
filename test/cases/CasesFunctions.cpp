#include "CasesFunctions.h"
#include <regex>
#include <nlohmann/json.hpp>

namespace maxwell{

    using namespace mfem;
    using json = nlohmann::json;



json parseJSONfile(const std::string& json_file)
{
	std::ifstream test_file(json_file);
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
    std::ifstream in(mesh_path);
    Mesh res(in, 0, 0, 0);
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
    FieldFactor[E][X] = 0.0; 
    FieldFactor[E][Y] = 0.0; 
    FieldFactor[E][Z] = 0.0; 
    FieldFactor[H][X] = 0.0; 
    FieldFactor[H][Y] = 0.0; 
    FieldFactor[H][Z] = 0.0; 
}

void ExcitationCoeffs::loadInitialPolarizationValues(const FieldType& ft, const Vector& pol)
{
    switch(ft){
        case E:
            if(pol[X] != 0.0){ FieldFactor[E][X] = pol[X];} else {FieldFactor[E][X] = 0.0;}
            if(pol[Y] != 0.0){ FieldFactor[E][Y] = pol[Y];} else {FieldFactor[E][Y] = 0.0;}
            if(pol[Z] != 0.0){ FieldFactor[E][Z] = pol[Z];} else {FieldFactor[E][Z] = 0.0;}
            break;
        case H:
            if(pol[X] != 0.0){ FieldFactor[H][X] = pol[X];} else {FieldFactor[H][X] = 0.0;}
            if(pol[Y] != 0.0){ FieldFactor[H][Y] = pol[Y];} else {FieldFactor[H][Y] = 0.0;}
            if(pol[Z] != 0.0){ FieldFactor[H][Z] = pol[Z];} else {FieldFactor[H][Z] = 0.0;}
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
                if (pol[X] != 0.0)     { FieldFactor[H][Z] = 1.0/std::sqrt(2.0); FieldFactor[H][Y] = -1.0/std::sqrt(2.0);}
                else if (pol[Y] != 0.0){ FieldFactor[H][X] = 1.0/std::sqrt(2.0); FieldFactor[H][Z] = -1.0/std::sqrt(2.0);}
                else                   { FieldFactor[H][Y] = 1.0/std::sqrt(2.0); FieldFactor[H][X] = -1.0/std::sqrt(2.0);}
                break;
            case H:
                if (pol[X] != 0.0)     { FieldFactor[E][Z] = -1.0/std::sqrt(2.0); FieldFactor[E][Y] =  1.0/std::sqrt(2.0);}
                else if (pol[Y] != 0.0){ FieldFactor[E][X] = -1.0/std::sqrt(2.0); FieldFactor[E][Z] =  1.0/std::sqrt(2.0);}
                else                   { FieldFactor[E][Y] = -1.0/std::sqrt(2.0); FieldFactor[E][X] =  1.0/std::sqrt(2.0);}
                break;
        }
    }
    else if (numExcitedDirections == 3){
        double normfactor = 0.0;
        switch(ft){
            case E:
                if (pol[X] != 0.0){ FieldFactor[H][Z] += pol[X]; FieldFactor[H][Y] -= pol[X];}
                if (pol[Y] != 0.0){ FieldFactor[H][X] += pol[Y]; FieldFactor[H][Z] -= pol[Y];}
                if (pol[Z] != 0.0){ FieldFactor[H][Y] += pol[Z]; FieldFactor[H][X] -= pol[Z];}
                normaliseExcitedVector(FieldFactor[H][X], FieldFactor[H][Y], FieldFactor[H][Z]);
                break;
            case H:
                if (pol[X] != 0.0){ FieldFactor[E][Z] -= pol[X]; FieldFactor[E][Y] += pol[X];}
                if (pol[Y] != 0.0){ FieldFactor[E][X] -= pol[Y]; FieldFactor[E][Z] += pol[Y];}
                if (pol[Z] != 0.0){ FieldFactor[E][Y] -= pol[Z]; FieldFactor[E][X] += pol[Z];}
                normaliseExcitedVector(FieldFactor[E][X], FieldFactor[E][Y], FieldFactor[E][Z]);
                break;
        }
    }
}

int extractRankNumberForMesh(const std::string& folder_name) {
    static std::regex rank_regex(R"(mesh_rank(\d+)(?:\.\d+)?)");
    std::smatch match;
    if (std::regex_match(folder_name, match, rank_regex)) {
        return std::stoi(match[1].str());
    }
    throw std::runtime_error("Invalid rank file name: " + folder_name);
}

int extractRankNumberForFolder(const std::string& folder_name) {
    static std::regex rank_regex(R"(rank_(\d+))");
    std::smatch match;
    if (std::regex_match(folder_name, match, rank_regex)) {
        return std::stoi(match[1].str());
    }
    throw std::runtime_error("Invalid rank file name: " + folder_name);
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

void RMSDataCalculator::loadNodepos(const std::string& data_path)
{
    std::filesystem::path probes_dir = std::filesystem::path(data_path) / "DomainSnapshopProbes";
    if (!std::filesystem::exists(probes_dir) || !std::filesystem::is_directory(probes_dir)) {
        throw std::runtime_error("No DomainSnapshopProbes were generated in this case. Rerun case with one.");
    }

    std::vector<std::filesystem::directory_entry> entries;

    for (const auto& entry : std::filesystem::directory_iterator(probes_dir)) {
        if (!entry.is_directory()) continue;

        std::string folder_name = entry.path().filename().string();
        if (folder_name == "meshes") continue;

        if (folder_name.rfind("rank_", 0) == 0) {
            entries.push_back(entry);
        }
    }
    
    std::sort(entries.begin(), entries.end(),
              [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b) {
                  return a.path().filename().string() < b.path().filename().string();
              });

    for (const auto& entry : entries) {
        std::string folder_name = entry.path().filename().string();
        int rank_no = extractRankNumberForFolder(folder_name);

        std::filesystem::path ex_file = entry.path() / "cycle_000000" / "Ex.gf";
        if (!std::filesystem::exists(ex_file)) {
            std::cerr << "Warning: Ex.gf missing in " << ex_file.parent_path() << "\n";
            continue;
        }

        if (rank_no >= (int)meshes_.size()) {
            throw std::runtime_error("Mesh not available for rank " + std::to_string(rank_no));
        }

        nodepos_[rank_no] =
            buildDoFPositions(loadGridFunctionFromFile(ex_file.string(), meshes_[rank_no]));
    }
}

void RMSDataCalculator::loadMeshes(const std::string& data_path)
{
    std::filesystem::path meshPath = data_path + "/DomainSnapshopProbes/meshes";
    if (!std::filesystem::exists(meshPath) || !std::filesystem::is_directory(meshPath)) {
        throw std::runtime_error("No DomainSnapshopProbes were generated in this case. Rerun case with one.");
    }

    std::vector<std::filesystem::directory_entry> entries;

    for (const auto& entry : std::filesystem::directory_iterator(meshPath)) {
        if (!entry.is_regular_file()) continue;

        const auto& fname = entry.path().filename().string();
        if (fname.rfind("mesh_rank", 0) == 0) {
            entries.push_back(entry);
        }
    }

    std::sort(entries.begin(), entries.end(),
              [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b) {
                  return extractRankNumberForMesh(a.path().filename().string()) <
                         extractRankNumberForMesh(b.path().filename().string());
              });

    int rank_no = 0;
    for (const auto& entry : entries) {
        meshes_[rank_no] = loadMeshFromFile(entry.path().string());
        rank_no++;
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
        return std::make_unique<TM2DSinusoidalMode>(case_data["sources"][0]["magnitude"]["modes"], box_size);
    }
    else{
        throw std::runtime_error("Currently unsupported FunctionType in buildFunctionByType.");
    }
}

std::unique_ptr<TimeFunction> buildDxFunctionByType(const json& case_data)
{
    if(getFunctionTypeFromJson(case_data) == FunctionType::Resonant){
        std::vector<int> modes = case_data["sources"][0]["magnitude"]["modes"];
        std::vector<double> box_size(modes.size()); //Assuming box size equal to 1.0 x 1.0 (x1.0 if 3D)
        for (auto v = 0; v < box_size.size(); v++){
            box_size[v] = 1.0;
        }
        return std::make_unique<TM2D_Dx_SinusoidalMode>(case_data["sources"][0]["magnitude"]["modes"], box_size);
    }
    else{
        throw std::runtime_error("Currently unsupported FunctionType in buildFunctionByType.");
    }
}

std::unique_ptr<TimeFunction> buildDyFunctionByType(const json& case_data)
{
    if(getFunctionTypeFromJson(case_data) == FunctionType::Resonant){
        std::vector<int> modes = case_data["sources"][0]["magnitude"]["modes"];
        std::vector<double> box_size(modes.size()); //Assuming box size equal to 1.0 x 1.0 (x1.0 if 3D)
        for (auto v = 0; v < box_size.size(); v++){
            box_size[v] = 1.0;
        }
        return std::make_unique<TM2D_Dy_SinusoidalMode>(case_data["sources"][0]["magnitude"]["modes"], box_size);
    }
    else{
        throw std::runtime_error("Currently unsupported FunctionType in buildFunctionByType.");
    }
}

void RMSDataCalculator::initFunctions(const std::string& json_file)
{
    auto case_data = parseJSONfile(json_file);
    function_ = buildFunctionByType(case_data);
    dx_function_ = buildDxFunctionByType(case_data);
    dy_function_ = buildDyFunctionByType(case_data);
}

GridFunction getGridFunction(const std::string& grid_path, Mesh& mesh)
{
    std::ifstream in(grid_path);
    GridFunction res(&mesh, in);
	return res;
}

double getTime(const std::filesystem::path& time_path)
{
    double res = 0.0;
    std::ifstream in(time_path);
    in >> res;
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

FieldScaleFactor getInitialExcitationCoefficients(const std::string& json_file)
{
    auto case_data = parseJSONfile(json_file);
    auto ft = assignFieldType(case_data["sources"][0]["field_type"]);
    auto pol = assemble3DVector(case_data["sources"][0]["polarization"]);
    pol /= pol.Norml2();

    FieldScaleFactor res;
    res[E][X] = 0.0; res[E][Y] = 0.0; res[E][Z] = 0.0;
    res[H][X] = 0.0; res[H][Y] = 0.0; res[H][Z] = 0.0;
    switch(ft){
        case E:
        res[E][X] = pol[0]; res[E][Y] = pol[1]; res[E][Z] = pol[2];
        break;
        case H:
        res[H][X] = pol[0]; res[H][Y] = pol[1]; res[H][Z] = pol[2];
        break;
    }
    return res;
}

std::string getFieldDirString(const FieldType& ft, const Direction& dir)
{
    switch(ft){
        case E:
            switch(dir){
                case X:
                return "/Ex.gf";
                case Y:
                return "/Ey.gf";
                case Z:
                return "/Ez.gf";
                default:
                throw std::runtime_error("Input direction in getFieldDirString is not valid.");
            }
        break;
        case H:
            switch(dir){
                case X:
                return "/Hx.gf";
                case Y:
                return "/Hy.gf";
                case Z:
                return "/Hz.gf";
                default:
                throw std::runtime_error("Input direction in getFieldDirString is not valid.");
            }
        break;
        default:
        throw std::runtime_error("Input FieldType in getFieldDirString is not valid.");
    }
}

double computeAverageElementSize(const std::string& data_path)
{

    std::filesystem::path sim_path = std::filesystem::path(data_path) / "SimulationStats";
    std::vector<double> element_sizes;

    for (const auto& entry : std::filesystem::directory_iterator(sim_path))
    {
        const auto& path = entry.path();
        if (path.filename().string().find("statistics_rank") == 0 && path.extension() == ".dat")
        {
            std::ifstream file(path);
            if (!file.is_open())
            {
                throw std::runtime_error("Failed to open file: " + entry.path().string() + "\n");
            }

            std::string line;
            while (std::getline(file, line))
            {
                const std::string key = "Average Element Size in Mesh:";
                auto pos = line.find(key);
                if (pos != std::string::npos)
                {
                    std::string value_str = line.substr(pos + key.size());
                    double value = std::stod(value_str);
                    element_sizes.push_back(value);
                    break;
                }
            }
        }
    }

    if (element_sizes.empty())
    {
        throw std::runtime_error("No element sizes found, check simulationstats output folder.");
    }

    double sum = 0.0;
    for (double val : element_sizes)
        sum += val;

    return sum / element_sizes.size();
}

RMSDataCalculator::RMSDataCalculator(const std::string& data_path, const std::string& json_file)
{
    loadMeshes(data_path);
    loadNodepos(data_path);
    initFunctions(json_file);
    
    ExcitationCoeffs excCoeff(json_file);
    FieldScaleFactor initialFactor = getInitialExcitationCoefficients(json_file);
    
    double rms = 0.0;
    int ndofs = getGlobalNdofs(nodepos_);
    int timeCount = 0;

    std::unique_ptr<FieldScaleFactor> scaleFactor;

    for (auto r = 0; r < meshes_.size(); r++)
    {
        std::string rank_string("rank_" + std::to_string(r));
        std::filesystem::path rankPath = data_path + "/DomainSnapshopProbes/" + rank_string;

        std::vector<std::filesystem::directory_entry> entries;
        for (auto& entry : std::filesystem::directory_iterator(rankPath)) {
            if (entry.is_directory()) { // optional: only subfolders
                entries.push_back(entry);
            }
        }

        std::sort(entries.begin(), entries.end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b){
            return a.path().filename() < b.path().filename();
        });

        for (auto& entry : entries) {

            std::filesystem::path cyclePath = entry.path();

            std::string timePath = cyclePath.string() + "/time.txt";
            double time = getTime(timePath);

            Vector analytic(nodepos_[r].size()); 

            for (auto v = 0; v < analytic.Size(); v++){
                analytic[v] = function_->eval(nodepos_[r][v], time);
            }

            time == 0.0 ? scaleFactor = std::make_unique<FieldScaleFactor>(initialFactor) : scaleFactor = std::make_unique<FieldScaleFactor>(excCoeff.FieldFactor);

            for (auto f : {E, H}){
				for (auto d : {X, Y, Z}){
                    std::string gfPath = cyclePath.string() + getFieldDirString(f, d);
                    auto gf = getGridFunction(gfPath, meshes_[r]);
                    double rootFactor = 0.0;
                    for (auto v = 0; v < analytic.Size(); v++){
                        if(f == E && d == Z){
                            rootFactor += std::pow(gf[v] - analytic[v] * scaleFactor->at(f)[d], 2.0);
                        }
                    }
                    rms += rootFactor;
                }
            }
            timeCount++;
        }
    }

    const auto& numRanks = meshes_.size();
    timeCount /= numRanks;

    double final_rms = std::sqrt(rms / (timeCount * ndofs));

    std::filesystem::path export_path = data_path + "/SimulationStats/AnalyticRMS.txt";
    if (!std::filesystem::exists(export_path.parent_path())) {
        if (!std::filesystem::create_directories(export_path.parent_path())) {
            throw std::runtime_error("Failed to create directory: " + export_path.parent_path().string());
        }
    }

    std::ofstream out_file(export_path);
    if (!out_file) {
        throw std::runtime_error("Error opening file for writing: " + export_path.string());
    }
    out_file << "RMS: " << final_rms << "\n";
    out_file << "Avg Element Size: " << computeAverageElementSize(data_path) << "\n";
    out_file.close();   

}


}