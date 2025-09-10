#include "CasesFunctions.h"

namespace maxwell{

    using namespace mfem;

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

int getRankAmount(const std::string& data_path)
{
    std::filesystem::path meshPath = std::filesystem::path(data_path) / "/DomainSnapshopProbes/meshes";
    if (!std::filesystem::exists(meshPath) || !std::filesystem::is_directory(meshPath)) {
        throw std::runtime_error("No DomainSnapshopProbes were generated in this case. Rerun case with one.");
    }

    int count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(meshPath)) {
        if (std::filesystem::is_regular_file(entry.status())) {
            count++;
        }
    }
    return count;
}

void L2SimDataCalculator::loadMeshes(const std::string& data_path)
{
    std::filesystem::path meshPath = std::filesystem::path(data_path) / "/DomainSnapshopProbes/meshes";
    if (!std::filesystem::exists(meshPath) || !std::filesystem::is_directory(meshPath)) {
        throw std::runtime_error("No DomainSnapshopProbes were generated in this case. Rerun case with one.");
    }

    std::vector<std::filesystem::path> meshFiles;

    for (const auto& entry : std::filesystem::directory_iterator(meshPath)) {
        if (entry.is_regular_file()) {
            const auto& fname = entry.path().filename().string();
            if (fname.rfind("mesh_rank", 0) == 0) {
                meshFiles.push_back(entry.path());
            }
        }
    }

    if (meshFiles.empty()) {
        throw std::runtime_error("No mesh files found in: " + meshPath.string());
    }

    for (const auto& file : meshFiles) {
        meshes_.push_back(loadMeshFromFile(file.string()));
    }
}

L2SimDataCalculator::L2SimDataCalculator(const std::string& data_path, const FunctionType function_type)
{
    loadMeshes(data_path);
}



}