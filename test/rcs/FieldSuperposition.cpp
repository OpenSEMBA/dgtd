#include "rcs/FieldSuperposition.h"

namespace maxwell{

CaseInfo::CaseInfo(const std::string& d, const std::string& j) :
data_path{d}, json_path{j} 
{};

void FieldSuperposition::loadRankPaths()
{
    c1_.rank_paths = findRankFolders(c1_.data_path);
    c2_.rank_paths = findRankFolders(c2_.data_path);
    c_compare_.rank_paths = findRankFolders(c_compare_.data_path);
    MFEM_ASSERT(c1_.rank_paths.size() == c2_.rank_paths.size() == c_compare_.rank_paths.size(), "Cases have different amounts of ranks.");
}

Mesh FieldSuperposition::loadCaseMesh(const size_t& rank)
{
    auto m1 = Mesh::LoadFromFile(c1_.rank_paths[rank], 1, 0, false);
    auto m2 = Mesh::LoadFromFile(c2_.rank_paths[rank], 1, 0, false);
    auto mc = Mesh::LoadFromFile(c_compare_.rank_paths[rank], 1, 0, false);
    MFEM_ASSERT(m1.GetNE() == m2.GetNE() == mc.GetNE(), "Meshes for different cases have mismatching number of elements.");
    MFEM_ASSERT(m1.GetNV() == m2.GetNV() == mc.GetNV(), "Meshes for different cases have mismatching number of vertices.");
    MFEM_ASSERT(m1.GetNBE() == m2.GetNBE() == mc.GetNBE(), "Meshes for different cases have mismatching number of boundary elements.");

    return m1;
}

FieldSuperposition::FieldSuperposition(const CaseInfo& c1, const CaseInfo& c2, const CaseInfo& c_compare) :
c1_(c1), c2_(c2), c_compare_(c_compare) 
{
    
    for (auto r = 0; r < c1_.rank_paths.size(); r++){
        auto mesh = loadCaseMesh(r);


    }




}






}