#include "rcs/FieldSuperposition.h"
#include "driver/driver.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <cmath>

namespace maxwell {

using namespace mfem;

struct AxisK {
    int axis;
    double k;
};

std::unique_ptr<GridFunction> initGridFunction(const GridFunction& ref)
{
    auto dgfec = dynamic_cast<const DG_FECollection*>(ref.FESpace()->FEColl());
    auto fes = FiniteElementSpace(ref.FESpace()->GetMesh(), dgfec);
    auto res = std::make_unique<GridFunction>(&fes);
    *res.get() = 0.0;
    return res;
}

std::vector<std::string> findRankFoldersOrdered(const std::string& root)
{
    std::vector<std::pair<int,std::string>> tmp;

    if (!std::filesystem::exists(root) || !std::filesystem::is_directory(root)) {
        throw std::runtime_error("data_path is not a directory: " + root);
    }

    std::regex pat(R"(rank_([0-9]+))");

    for (const auto& de : std::filesystem::directory_iterator(root)) {
        if (!de.is_directory()) continue;
        const std::string name = de.path().filename().string();
        std::smatch m;
        if (std::regex_match(name, m, pat)) {
            int idx = std::stoi(m[1].str());
            tmp.emplace_back(idx, de.path().string());
        }
    }

    std::sort(tmp.begin(), tmp.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    if (tmp.empty() || tmp.front().first != 0) {
        throw std::runtime_error("missing required rank_0 folder in " + root);
    }

    std::vector<std::string> out;
    out.reserve(tmp.size());
    for (auto& p : tmp) out.push_back(std::move(p.second));
    return out;
}

CaseInfo::CaseInfo(const std::string& d, const std::string& j)
    : data_path{d}, json_path{j} {}

double readTimeFile(const std::filesystem::path& snapshot_dir)
{
    const std::filesystem::path time_path = snapshot_dir / "time.txt";
    std::ifstream in(time_path);
    if (!in) { throw std::runtime_error("Missing time.txt in " + snapshot_dir.string()); }
    double t = 0.0;
    in >> t;
    if (!in) { throw std::runtime_error("Invalid time.txt in " + snapshot_dir.string()); }
    return t;
}

bool checkHasAllFieldFiles(const std::filesystem::path& snapshot_dir)
{
    const std::vector<std::string> names = {"Ex.gf","Ey.gf","Ez.gf","Hx.gf","Hy.gf","Hz.gf","time.txt"};
    for (const auto& n : names) {
        if (!std::filesystem::exists(snapshot_dir / n)) { return false; }
    }
    return true;
}

std::unique_ptr<GridFunction> loadGridFunction(const std::filesystem::path& file, Mesh& mesh)
{
    std::ifstream in(file, std::ios::binary);
    if (!in) { throw std::runtime_error("Cannot open " + file.string()); }
    auto res = std::make_unique<GridFunction>(&mesh, in);
    return res;
}

std::vector<std::pair<double, std::filesystem::path>> getAndReorderSnapshots(const std::filesystem::path& base_dir)
{
    std::vector<std::pair<double, std::filesystem::path>> time_dirs;
    if (!std::filesystem::exists(base_dir) || !std::filesystem::is_directory(base_dir)) {
        throw std::runtime_error("Base directory not found: " + base_dir.string());
    }

    for (const auto& entry : std::filesystem::directory_iterator(base_dir)) {
        if (!entry.is_directory()) { continue; }
        const std::filesystem::path snap = entry.path();
        if (!checkHasAllFieldFiles(snap)) { continue; }
        const auto time = readTimeFile(snap);
        time_dirs.emplace_back(time, snap);
    }

    std::sort(time_dirs.begin(), time_dirs.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });
    return time_dirs;
}

TimeFields loadSnapshotFields(const std::filesystem::path& snapshot_dir, Mesh& mesh)
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

TimeToFields buildTimeToFields(const std::string& rank_root, Mesh& mesh)
{
    TimeToFields t2f;
    const auto ordered = getAndReorderSnapshots(rank_root);
    for (const auto& [t, dir] : ordered) {
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
    MFEM_ASSERT(c1_.rank_paths.size() == c2_.rank_paths.size() &&
                c2_.rank_paths.size() == c_compare_.rank_paths.size(),
                "Cases have different amounts of ranks.");
}

std::string findMeshFileForRank(const std::string& data_path, size_t rank)
{
    const std::filesystem::path meshes_dir = std::filesystem::path(data_path) / "meshes";
    if (!std::filesystem::exists(meshes_dir) || !std::filesystem::is_directory(meshes_dir)) {
        throw std::runtime_error("meshes directory not found: " + meshes_dir.string());
    }
    const std::string prefix_with_us = "mesh_rank_" + std::to_string(rank) + ".";
    const std::string prefix_no_us   = "mesh_rank"  + std::to_string(rank) + ".";
    std::vector<std::filesystem::path> matches;
    for (const auto& de : std::filesystem::directory_iterator(meshes_dir)) {
        if (!de.is_regular_file()) continue;
        const std::string fname = de.path().filename().string();
        if (fname.rfind(prefix_with_us, 0) == 0 || fname.rfind(prefix_no_us, 0) == 0) {
            matches.push_back(de.path());
        }
    }
    if (matches.empty()) {
        throw std::runtime_error("no mesh file for rank " + std::to_string(rank) + " under " + meshes_dir.string());
    }
    std::sort(matches.begin(), matches.end());
    return matches.back().string();
}

Mesh FieldSuperposition::loadCaseMesh(const size_t& rank)
{
    const std::string m1p = findMeshFileForRank(c1_.data_path, rank);
    const std::string m2p = findMeshFileForRank(c2_.data_path, rank);
    const std::string mcp = findMeshFileForRank(c_compare_.data_path, rank);

    auto m1 = Mesh::LoadFromFile(m1p.c_str(), 1, 0, false);
    auto m2 = Mesh::LoadFromFile(m2p.c_str(), 1, 0, false);
    auto mc = Mesh::LoadFromFile(mcp.c_str(), 1, 0, false);
    MFEM_ASSERT(m1.GetNE() == m2.GetNE() && m2.GetNE() == mc.GetNE(),  "Meshes: NE mismatch");
    MFEM_ASSERT(m1.GetNV() == m2.GetNV() && m2.GetNV() == mc.GetNV(),  "Meshes: NV mismatch");
    MFEM_ASSERT(m1.GetNBE()== m2.GetNBE() && m2.GetNBE()== mc.GetNBE(), "Meshes: NBE mismatch");
    return m1;
}

const GridFunction& pickComponent(const TimeFields& tf, const char* name)
{
    if      (std::strcmp(name,"Ex")==0) return *tf.Ex;
    else if (std::strcmp(name,"Ey")==0) return *tf.Ey;
    else if (std::strcmp(name,"Ez")==0) return *tf.Ez;
    else if (std::strcmp(name,"Hx")==0) return *tf.Hx;
    else if (std::strcmp(name,"Hy")==0) return *tf.Hy;
    else if (std::strcmp(name,"Hz")==0) return *tf.Hz;
    else { throw std::runtime_error("Unknown component passed as name in pickComponent."); }
}

ComplexField addFrequencyComponentFromCache(const TimeToFields& t2f,
                                            const char* component,
                                            double freq_hz)
{
    MFEM_VERIFY(!t2f.empty(), "empty time series");

    const GridFunction& ref = pickComponent(t2f.begin()->second, component);

    ComplexField out;
    out.real = initGridFunction(ref);
    out.imag = initGridFunction(ref);

    const int nd = ref.Size();

    for (const auto& kv : t2f) {
        const double t = kv.first;
        const TimeFields& tf = kv.second;
        const GridFunction& gf = pickComponent(tf, component);

        MFEM_VERIFY(gf.Size() == nd, "size mismatch across time snapshots");

        const double arg = 2.0 * M_PI * freq_hz * t;
        const double cos_factor = std::cos(arg);
        const double sin_factor = std::sin(arg);

        for (int i = 0; i < nd; ++i) {
            const double x = gf[i];
            (*out.real)[i] += x * cos_factor;
            (*out.imag)[i] += -x * sin_factor;
        }
    }

    return out;
}

FrequencyFields addAllFieldsAtFrequencyFromCache(const TimeToFields& t2f, double freq_hz, FiniteElementSpace& fes)
{
    FrequencyFields F;
    F.Ex = addFrequencyComponentFromCache(t2f, "Ex", freq_hz);
    F.Ey = addFrequencyComponentFromCache(t2f, "Ey", freq_hz);
    F.Ez = addFrequencyComponentFromCache(t2f, "Ez", freq_hz);
    F.Hx = addFrequencyComponentFromCache(t2f, "Hx", freq_hz);
    F.Hy = addFrequencyComponentFromCache(t2f, "Hy", freq_hz);
    F.Hz = addFrequencyComponentFromCache(t2f, "Hz", freq_hz);
    F.Ex.real.get()->SetSpace(&fes);
    F.Ey.real.get()->SetSpace(&fes);
    F.Ez.real.get()->SetSpace(&fes);
    F.Hx.real.get()->SetSpace(&fes);
    F.Hy.real.get()->SetSpace(&fes);
    F.Hz.real.get()->SetSpace(&fes);
    F.Ex.imag.get()->SetSpace(&fes);
    F.Ey.imag.get()->SetSpace(&fes);
    F.Ez.imag.get()->SetSpace(&fes);
    F.Hx.imag.get()->SetSpace(&fes);
    F.Hy.imag.get()->SetSpace(&fes);
    F.Hz.imag.get()->SetSpace(&fes);
    return F;
}

Vector loadPropagationK(const std::string& json_path)
{
    const auto case_data = driver::parseJSONfile(json_path);
    return driver::assemble3DVector(case_data["sources"][0]["propagation"]);
}

AxisK selectAxisAndK(const Vector& kvec)
{
    int axis = 0;
    double mags[3] = { std::abs(kvec[0]), std::abs(kvec[1]), std::abs(kvec[2]) };
    if (mags[1] > mags[axis]) axis = 1;
    if (mags[2] > mags[axis]) axis = 2;
    return { axis, kvec[axis] };
}

void applyPhaseInPlace(ComplexField& field, const GridFunction& nodes, int axis, double k)
{
    const int ndofs = field.real->Size();
    MFEM_ASSERT(field.imag->Size() == ndofs, "ComplexField size mismatch");
    MFEM_ASSERT(nodes.VectorDim() >= axis+1, "Nodes missing requested axis");
    MFEM_ASSERT(nodes.Size() >= (axis+1)*ndofs, "Nodes size inconsistent");

    const int offset = axis * ndofs;

    for (int i = 0; i < ndofs; ++i) {
        const double coord = nodes[offset + i];
        const double phase = k * coord;
        const double cos_factor = std::cos(phase);
        const double sin_factor = std::sin(phase);
        const double real_part = (*field.real)[i];
        const double imag_part = (*field.imag)[i];
        (*field.real)[i] =  real_part*cos_factor - imag_part*sin_factor;
        (*field.imag)[i] =  real_part*sin_factor + imag_part*cos_factor;
    }
}

void applyPhaseInPlace(ComplexField& field, const std::vector<mfem::Vector>& dof_pos, int axis, double k)
{
    const int ndofs = field.real->Size();
    MFEM_ASSERT(field.imag->Size() == ndofs, "ComplexField size mismatch");
    MFEM_ASSERT(static_cast<int>(dof_pos.size()) == ndofs, "DOF positions size mismatch");

    for (int i = 0; i < ndofs; ++i) {
        const double coord = dof_pos[i](axis);
        const double phase = k * coord;
        const double cos_factor = std::cos(phase);
        const double sin_factor = std::sin(phase);
        const double real_part = (*field.real)[i];
        const double imag_part = (*field.imag)[i];
        (*field.real)[i] =  real_part*cos_factor - imag_part*sin_factor;
        (*field.imag)[i] =  real_part*sin_factor + imag_part*cos_factor;
    }
}

void applyPhaseToFrequencyFields(FrequencyFields& freqFields, const GridFunction& nodes, int axis, double k)
{
    applyPhaseInPlace(freqFields.Ex, nodes, axis, k);
    applyPhaseInPlace(freqFields.Ey, nodes, axis, k);
    applyPhaseInPlace(freqFields.Ez, nodes, axis, k);
    applyPhaseInPlace(freqFields.Hx, nodes, axis, k);
    applyPhaseInPlace(freqFields.Hy, nodes, axis, k);
    applyPhaseInPlace(freqFields.Hz, nodes, axis, k);
}

void applyPhaseToFrequencyFields(FrequencyFields& freqFields, const std::vector<mfem::Vector>& dof_pos, int axis, double k)
{
    applyPhaseInPlace(freqFields.Ex, dof_pos, axis, k);
    applyPhaseInPlace(freqFields.Ey, dof_pos, axis, k);
    applyPhaseInPlace(freqFields.Ez, dof_pos, axis, k);
    applyPhaseInPlace(freqFields.Hx, dof_pos, axis, k);
    applyPhaseInPlace(freqFields.Hy, dof_pos, axis, k);
    applyPhaseInPlace(freqFields.Hz, dof_pos, axis, k);
}

std::vector<Vector> buildDofPositions(const FiniteElementSpace& fes)
{
    auto fec = dynamic_cast<const DG_FECollection*>(fes.FEColl());
    MFEM_VERIFY(fec != nullptr, "Expected DG_FECollection for DOF positions");
    FiniteElementSpace vdimfes(fes.GetMesh(), fec, 3);
    GridFunction nodes(&vdimfes);
    fes.GetMesh()->GetNodes(nodes);
    const int dirSize = nodes.Size() / 3;
    std::vector<mfem::Vector> res;
    res.reserve(dirSize);
    for (int i = 0; i < dirSize; ++i) {
        mfem::Vector p(3);
        p[0] = nodes[i];
        p[1] = nodes[i + dirSize];
        p[2] = nodes[i + 2*dirSize];
        res.emplace_back(std::move(p));
    }
    return res;
}

void exportFrequencyGridFunctions(const FrequencyFields& ff, const CaseInfo& ci, size_t rank, const std::string& stage)
{
    std::filesystem::path base(ci.data_path);
    std::filesystem::path parent = base;
    if (parent.filename().empty()) {
        parent = base.parent_path();
    }
    parent = parent.parent_path();
    std::filesystem::path out = parent / "FieldSuperposition" / ("rank_" + std::to_string(rank)) / stage;
    std::filesystem::create_directories(out);

    Mesh* mesh = ff.Ex.real->FESpace()->GetMesh();

    std::string paraview_case_name = "fsp_" + parent.filename().string();
    ParaViewDataCollection dc(paraview_case_name, mesh);
    dc.SetPrefixPath(out.string());
    dc.SetHighOrderOutput(true);

    dc.RegisterField("Ex_real", const_cast<GridFunction*>(ff.Ex.real.get()));
    dc.RegisterField("Ex_imag", const_cast<GridFunction*>(ff.Ex.imag.get()));
    dc.RegisterField("Ey_real", const_cast<GridFunction*>(ff.Ey.real.get()));
    dc.RegisterField("Ey_imag", const_cast<GridFunction*>(ff.Ey.imag.get()));
    dc.RegisterField("Ez_real", const_cast<GridFunction*>(ff.Ez.real.get()));
    dc.RegisterField("Ez_imag", const_cast<GridFunction*>(ff.Ez.imag.get()));
    dc.RegisterField("Hx_real", const_cast<GridFunction*>(ff.Hx.real.get()));
    dc.RegisterField("Hx_imag", const_cast<GridFunction*>(ff.Hx.imag.get()));
    dc.RegisterField("Hy_real", const_cast<GridFunction*>(ff.Hy.real.get()));
    dc.RegisterField("Hy_imag", const_cast<GridFunction*>(ff.Hy.imag.get()));
    dc.RegisterField("Hz_real", const_cast<GridFunction*>(ff.Hz.real.get()));
    dc.RegisterField("Hz_imag", const_cast<GridFunction*>(ff.Hz.imag.get()));

    dc.Save();
}

FieldSuperposition::FieldSuperposition(const CaseInfo& c1, const CaseInfo& c2, const CaseInfo& c_compare, const Frequency freq)
: c1_(c1), c2_(c2), c_compare_(c_compare)
{
    c1_.rank_paths = findRankFoldersOrdered(c1_.data_path);
    c2_.rank_paths = findRankFoldersOrdered(c2_.data_path);
    c_compare_.rank_paths = findRankFoldersOrdered(c_compare_.data_path);

    if (!(c1_.rank_paths.size() == c2_.rank_paths.size() && c2_.rank_paths.size() == c_compare_.rank_paths.size())) {
        throw std::runtime_error("Number of ranks for simulations to compare mismatches, verify exported simulation data.");
    }

    for (size_t r = 0; r < c1_.rank_paths.size(); ++r) {
        Mesh mesh = loadCaseMesh(r);

        TimeToFields t2f1 = buildTimeToFields((std::filesystem::path(c1_.data_path) / ("rank_" + std::to_string(r))).string(), mesh);
        TimeToFields t2f2 = buildTimeToFields((std::filesystem::path(c2_.data_path) / ("rank_" + std::to_string(r))).string(), mesh);
        TimeToFields t2fc = buildTimeToFields((std::filesystem::path(c_compare_.data_path) / ("rank_" + std::to_string(r))).string(), mesh);

        auto fes = t2f1.begin()->second.Ex->FESpace();

        FrequencyFields ff1 = addAllFieldsAtFrequencyFromCache(t2f1, freq, *fes);
        FrequencyFields ff2 = addAllFieldsAtFrequencyFromCache(t2f2, freq, *fes);
        FrequencyFields ff_compare = addAllFieldsAtFrequencyFromCache(t2fc, freq, *fes);

        exportFrequencyGridFunctions(ff1, c1_, r, "precoeff");
        exportFrequencyGridFunctions(ff2, c2_, r, "precoeff");
        exportFrequencyGridFunctions(ff_compare, c_compare_, r, "precoeff");

        const GridFunction& ref_gf = *t2f1.begin()->second.Ex;
        std::vector<mfem::Vector> dof_pos = buildDofPositions(*ref_gf.FESpace());

        const Vector k1 = loadPropagationK(c1_.json_path);
        const Vector k2 = loadPropagationK(c2_.json_path);
        const AxisK a1 = selectAxisAndK(k1);
        const AxisK a2 = selectAxisAndK(k2);

        applyPhaseToFrequencyFields(ff1, dof_pos, a2.axis, a2.k);
        applyPhaseToFrequencyFields(ff2, dof_pos, a1.axis, a1.k);

        exportFrequencyGridFunctions(ff1, c1_, r, "postcoeff");
        exportFrequencyGridFunctions(ff2, c2_, r, "postcoeff");
        exportFrequencyGridFunctions(ff_compare, c_compare_, r, "postcoeff");
    }
}

}
