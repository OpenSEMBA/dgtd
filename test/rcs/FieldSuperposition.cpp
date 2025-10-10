#include "rcs/FieldSuperposition.h"
#include "driver/driver.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <algorithm>

namespace maxwell{

using namespace mfem;

CaseInfo::CaseInfo(const std::string& d, const std::string& j) :
    data_path{d}, json_path{j} 
{};


using Frequency = double;
using FrequencyToFields = std::map<Frequency, FrequencyFields>;

double readTimeFile(const std::filesystem::path &snapshot_dir)
{
    const std::filesystem::path tf = snapshot_dir / "time.txt";
    std::ifstream in(tf);
    if (!in) { throw std::runtime_error("Missing time.txt in " + snapshot_dir.string()); }
    double t = 0.0;
    in >> t;
    if (!in) { throw std::runtime_error("Invalid time.txt in " + snapshot_dir.string()); }
    return t;
}

bool checkHasAllFieldFiles(const std::filesystem::path &snapshot_dir)
{
    static const char* names[] = {"Ex.gf","Ey.gf","Ez.gf","Hx.gf","Hy.gf","Hz.gf","time.txt"};
    for (const char* n : names) {
        if (!std::filesystem::exists(snapshot_dir / n)) { return false; }
    }
    return true;
}

std::unique_ptr<GridFunction> loadGridFunction(const std::filesystem::path &file, Mesh &mesh)
{
    std::ifstream in(file, std::ios::binary);
    if (!in) { 
        throw std::runtime_error("Cannot open " + file.string()); 
    }
    auto res = std::make_unique<GridFunction>(&mesh, in);
    return res;
}

std::vector<std::pair<double, std::filesystem::path>> getAndReorderSnapshots(const std::filesystem::path &base_dir)
{
    std::vector<std::pair<double, std::filesystem::path>> time_dirs;
    if (!std::filesystem::exists(base_dir) || !std::filesystem::is_directory(base_dir)) {
        throw std::runtime_error("Base directory not found: " + base_dir.string());
    }

    for (const auto &entry : std::filesystem::directory_iterator(base_dir)) {
        if (!entry.is_directory()) { continue; }
        const std::filesystem::path snap = entry.path();
        if (!checkHasAllFieldFiles(snap)) { continue; }
        const double t = readTimeFile(snap);
        time_dirs.emplace_back(t, snap);
    }

    std::sort(time_dirs.begin(), time_dirs.end(),
                [](const auto &a, const auto &b){ return a.first < b.first; });
    return time_dirs;
}

TimeFields loadSnapshotFields(const std::filesystem::path &snapshot_dir, Mesh &mesh)
{
    TimeFields tfs;
    tfs.Ex = loadGridFunction(snapshot_dir / "Ex.gf", mesh);
    tfs.Ey = loadGridFunction(snapshot_dir / "Ey.gf", mesh);
    tfs.Ez = loadGridFunction(snapshot_dir / "Ez.gf", mesh);
    tfs.Hx = loadGridFunction(snapshot_dir / "Hx.gf", mesh);
    tfs.Hy = loadGridFunction(snapshot_dir / "Hy.gf", mesh);
    tfs.Hz = loadGridFunction(snapshot_dir / "Hz.gf", mesh);
    return tfs;
}

TimeToFields buildTimeToFields(const std::string &rank_root, mfem::Mesh &mesh)
{
    TimeToFields t2f;

    const auto ordered = getAndReorderSnapshots(rank_root);
    for (const auto &[t, dir] : ordered) {
        TimeFields tfs = loadSnapshotFields(dir, mesh);
        t2f.emplace(t, std::move(tfs));
    }
    return t2f;
}

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


inline std::unique_ptr<GridFunction> initGridFunction(const GridFunction& ref)
{
    auto dgfec = dynamic_cast<const DG_FECollection*>(ref.FESpace()->FEColl());
    auto fes = FiniteElementSpace(ref.FESpace()->GetMesh(), dgfec);
    auto res = std::make_unique<GridFunction>(&fes);
    *res.get() = 0.0;
    return res;
}

inline ComplexField addFrequencyComponent(Mesh& mesh, const std::string& base_path, const char* filename, const double freq)
{

    std::vector<std::pair<double, std::filesystem::path>> snaps = getAndReorderSnapshots(base_path);
    MFEM_ASSERT(!snaps.empty(), "No snapshots found in " + base_path);

    ComplexField res;
    bool inited = false;

    for (const auto &[t, dir] : snaps) {
        auto gf_ptr = loadGridFunction(dir / filename, mesh);
        const GridFunction& gf = *gf_ptr;

        if (!inited) {
        res.re = initGridFunction(gf);
        res.im = initGridFunction(gf);
        inited = true;
        }

        MFEM_ASSERT(res.re->Size() == gf.Size() && res.im->Size() == gf.Size(), "Size mismatch during accumulation");

        const double arg = 2.0 * M_PI * freq * t;
        const double real_exp = std::cos(arg);
        const double imag_exp = std::sin(arg);

        for (int i = 0; i < gf.Size(); ++i) {
        res.re.get()[i] += gf[i] * real_exp;
        res.im.get()[i] -= gf[i] * imag_exp;
        }
    }

    //normalization by total nodes would go here if needed.

    return res;
}

inline FrequencyFields addAllFieldsAtFrequency(Mesh &mesh, const std::string &base_path, double freq)
{
    FrequencyFields res;
    res.Ex = addFrequencyComponent(mesh, base_path, "Ex.gf", freq);
    res.Ey = addFrequencyComponent(mesh, base_path, "Ey.gf", freq);
    res.Ez = addFrequencyComponent(mesh, base_path, "Ez.gf", freq);
    res.Hx = addFrequencyComponent(mesh, base_path, "Hx.gf", freq);
    res.Hy = addFrequencyComponent(mesh, base_path, "Hy.gf", freq);
    res.Hz = addFrequencyComponent(mesh, base_path, "Hz.gf", freq);
    return res;
}

FieldSuperposition::FieldSuperposition(const CaseInfo& c1, const CaseInfo& c2, const CaseInfo& c_compare, const Frequency freq) :
c1_(c1), c2_(c2), c_compare_(c_compare) 
{
    for (auto r = 0; r < c1_.rank_paths.size(); r++){
        auto mesh = loadCaseMesh(r);
        TimeToFields t2f1_ = buildTimeToFields(c1_.rank_paths[r], mesh);
        TimeToFields t2f2_ = buildTimeToFields(c2_.rank_paths[r], mesh);
        TimeToFields t2f_compare_ = buildTimeToFields(c_compare_.rank_paths[r], mesh);

        FrequencyFields ff1 = addAllFieldsAtFrequency(mesh, c1_.rank_paths[r], freq);
        FrequencyFields ff2 = addAllFieldsAtFrequency(mesh, c2_.rank_paths[r], freq);
        FrequencyFields ff_compare = addAllFieldsAtFrequency(mesh, c_compare_.rank_paths[r], freq);

    }

}






}