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

CaseInfo::CaseInfo(const std::string& d, const std::string& j)
    : data_path{d}, json_path{j} {}


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
    if (!in) { throw std::runtime_error("Cannot open " + file.string()); }
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

TimeToFields buildTimeToFields(const std::string &rank_root, Mesh &mesh)
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
    MFEM_ASSERT(c1_.rank_paths.size() == c2_.rank_paths.size() &&
                c2_.rank_paths.size() == c_compare_.rank_paths.size(),
                "Cases have different amounts of ranks.");
}

Mesh FieldSuperposition::loadCaseMesh(const size_t& rank)
{
    auto m1 = Mesh::LoadFromFile(c1_.rank_paths[rank], 1, 0, false);
    auto m2 = Mesh::LoadFromFile(c2_.rank_paths[rank], 1, 0, false);
    auto mc = Mesh::LoadFromFile(c_compare_.rank_paths[rank], 1, 0, false);
    MFEM_ASSERT(m1.GetNE() == m2.GetNE() && m2.GetNE() == mc.GetNE(),  "Meshes: NE mismatch");
    MFEM_ASSERT(m1.GetNV() == m2.GetNV() && m2.GetNV() == mc.GetNV(),  "Meshes: NV mismatch");
    MFEM_ASSERT(m1.GetNBE()== m2.GetNBE()&& m2.GetNBE()== mc.GetNBE(), "Meshes: NBE mismatch");
    return m1;
}

const GridFunction& pickComponent(const TimeFields &tf, const char *name)
{
    if      (std::strcmp(name,"Ex")==0) return *tf.Ex;
    else if (std::strcmp(name,"Ey")==0) return *tf.Ey;
    else if (std::strcmp(name,"Ez")==0) return *tf.Ez;
    else if (std::strcmp(name,"Hx")==0) return *tf.Hx;
    else if (std::strcmp(name,"Hy")==0) return *tf.Hy;
    else if (std::strcmp(name,"Hz")==0) return *tf.Hz;
    else{
        throw std::runtime_error("Unknown component passed as name in pickComponent.");
    }
}

ComplexField addFrequencyComponentFromCache(const TimeToFields &t2f,
                                                   const char *component,
                                                   double freq_hz)
{
    MFEM_VERIFY(!t2f.empty(), "empty time series");

    const GridFunction &ref = pickComponent(t2f.begin()->second, component);

    ComplexField out;
    out.real = initGridFunction(ref);
    out.imag = initGridFunction(ref);

    const int nd = ref.Size();

    for (const auto &kv : t2f) {
        const double t = kv.first;
        const TimeFields &tf = kv.second;
        const GridFunction &gf = pickComponent(tf, component);

        MFEM_VERIFY(gf.Size() == nd, "size mismatch across time snapshots");

        const double arg = 2.0 * M_PI * freq_hz * t;
        const double c = std::cos(arg);
        const double s = std::sin(arg);

        for (int i = 0; i < nd; ++i) {
            const double x = gf[i];
            (*out.real)[i] += x * c;
            (*out.imag)[i] += -x * s; // e^{-jωt}
        }
    }

    // optional normalization goes here.

    return out;
}

FrequencyFields addAllFieldsAtFrequencyFromCache(const TimeToFields &t2f, double freq_hz)
{
    FrequencyFields F;
    F.Ex = addFrequencyComponentFromCache(t2f, "Ex", freq_hz);
    F.Ey = addFrequencyComponentFromCache(t2f, "Ey", freq_hz);
    F.Ez = addFrequencyComponentFromCache(t2f, "Ez", freq_hz);
    F.Hx = addFrequencyComponentFromCache(t2f, "Hx", freq_hz);
    F.Hy = addFrequencyComponentFromCache(t2f, "Hy", freq_hz);
    F.Hz = addFrequencyComponentFromCache(t2f, "Hz", freq_hz);
    return F;
}

Vector loadPropagationK(const std::string &json_path)
{
    const auto case_data = driver::parseJSONfile(json_path);
    return driver::assemble3DVector(case_data["sources"][0]["propagation"]);
}

AxisK selectAxisAndK(const Vector &kvec)
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
        const double c = std::cos(phase);
        const double s = std::sin(phase);
        const double re0 = (*field.real)[i];
        const double im0 = (*field.imag)[i];
        (*field.real)[i] =  re0*c - im0*s;
        (*field.imag)[i] =  re0*s + im0*c;
    }
}

void applyPhaseToFrequencyFields(FrequencyFields &freqFields, const GridFunction &nodes, int axis, double k)
{
    applyPhaseInPlace(freqFields.Ex, nodes, axis, k);
    applyPhaseInPlace(freqFields.Ey, nodes, axis, k);
    applyPhaseInPlace(freqFields.Ez, nodes, axis, k);
    applyPhaseInPlace(freqFields.Hx, nodes, axis, k);
    applyPhaseInPlace(freqFields.Hy, nodes, axis, k);
    applyPhaseInPlace(freqFields.Hz, nodes, axis, k);
}

void crossApplyPhase(FrequencyFields &freq_1, FrequencyFields &freq_2, Mesh &mesh, const std::string &c1_json, const std::string &c2_json)
{
    mesh.EnsureNodes();
    MFEM_ASSERT(mesh.GetNodes(), "Mesh has no nodes");
    const GridFunction &nodes = *mesh.GetNodes();

    const Vector k1 = loadPropagationK(c1_json);
    const Vector k2 = loadPropagationK(c2_json);

    const AxisK a1 = selectAxisAndK(k1);
    const AxisK a2 = selectAxisAndK(k2);

    applyPhaseToFrequencyFields(freq_1, nodes, a2.axis, a2.k);
    applyPhaseToFrequencyFields(freq_2, nodes, a1.axis, a1.k);
}

FieldSuperposition::FieldSuperposition(const CaseInfo& c1, const CaseInfo& c2, const CaseInfo& c_compare, const Frequency freq)
: c1_(c1), c2_(c2), c_compare_(c_compare)
{
    for (size_t r = 0; r < c1_.rank_paths.size(); ++r) {
        Mesh mesh = loadCaseMesh(r);

        TimeToFields t2f1 = buildTimeToFields(c1_.rank_paths[r], mesh);
        TimeToFields t2f2 = buildTimeToFields(c2_.rank_paths[r], mesh);
        TimeToFields t2fc = buildTimeToFields(c_compare_.rank_paths[r], mesh);

        FrequencyFields ff1 = addAllFieldsAtFrequencyFromCache(t2f1, freq);
        FrequencyFields ff2 = addAllFieldsAtFrequencyFromCache(t2f2, freq);
        FrequencyFields ff_compare = addAllFieldsAtFrequencyFromCache(t2fc, freq);

        crossApplyPhase(ff1, ff2, mesh, c1_.json_path, c2_.json_path);

    }
}

}
