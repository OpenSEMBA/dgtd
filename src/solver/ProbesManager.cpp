#include "ProbesManager.h"
#include "SourcesManager.h"
#include "math/PhysicalConstants.h"
#include "general/text.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <fcntl.h>
#include <unistd.h>

namespace maxwell {

using namespace mfem;

namespace {

MPI_Comm getFESComm(const ParFiniteElementSpace& fes)
{
    return fes.GetParMesh()->GetComm();
}

bool isNodeRoot(MPI_Comm comm)
{
    MPI_Comm node_comm;
    int comm_rank = 0;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, comm_rank, MPI_INFO_NULL, &node_comm);
    
    int node_rank;
    MPI_Comm_rank(node_comm, &node_rank);
    MPI_Comm_free(&node_comm);
    
    return (node_rank == 0);
}

#ifdef SHOW_TIMER_INFORMATION
void formatProbeTimingCell(char* buf, std::size_t n, double ms, bool defined, int width)
{
    if (!defined) {
        std::snprintf(buf, n, "%*s", width, "off");
        return;
    }
    std::snprintf(buf, n, "%*.4f", width, ms);
}
#endif

}  // namespace

std::string getRunModeTag()
{
    std::string backend;
    if (mfem::Device::Allows(mfem::Backend::CUDA)){
        backend = "cuda-";
        backend.append(std::to_string(Mpi::WorldSize()));
        return backend;
    }
    else{
        if (Mpi::WorldSize() == 1){
            return "single-core";
        }
        else{
            backend = "mpi-";
            backend.append(std::to_string(Mpi::WorldSize()));
            return backend;
        }
    }
}

std::string getSimulationCaseExportPath(const std::string& caseName)
{
    return "exports/SimulationData/" + getRunModeTag() + "/" + caseName;
}

std::string getFieldPolString(const FieldType& ft, const Direction& d)
{
    switch(ft){
        case E:
        switch(d){
            case X: return "Ex";
            case Y: return "Ey";
            case Z: return "Ez";
            default: throw std::runtime_error("Incorrect direction in getFieldPolString.");
        }
        case H:
        switch(d){
            case X: return "Hx";
            case Y: return "Hy";
            case Z: return "Hz";
            default: throw std::runtime_error("Incorrect direction in getFieldPolString.");
        }
        default:
            throw std::runtime_error("Incorrect fieldtype in getFieldPolString.");
    }
}

ParaViewDataCollection ProbesManager::buildParaviewDataCollectionInfo(const ExporterProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{
    fes_.ExchangeFaceNbrData();
    fes_.GetParMesh()->ExchangeFaceNbrData();
    ParaViewDataCollection pd{ p.name, fes_.GetParMesh()};
    const MPI_Comm comm = getFESComm(fes_);
    
    std::string paraview_path = "exports/ParaView/" + getRunModeTag() + "/";
    
    if (isNodeRoot(comm)) {
        std::filesystem::create_directories(paraview_path);
    }
    MPI_Barrier(comm);

    pd.SetPrefixPath(paraview_path);

    pd.RegisterField("E", &fields.get(E));
    pd.RegisterField("H", &fields.get(H));
    
    bool highOrder = false;
    auto geomElemOrder = fes_.GetMesh()->GetElementTransformation(0)->Order();
    auto fecorder = fes_.FEColl()->GetOrder();
    geomElemOrder > 1 || fecorder > 1 ? highOrder = true : highOrder = false;
    pd.SetHighOrderOutput(highOrder);
    pd.SetLevelsOfDetail(std::max(geomElemOrder,fecorder));
    
    pd.SetDataFormat(VTKFormat::BINARY);

    return pd;
}

ProbesManager::ProbesManager(Probes pIn, mfem::ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>& fields, const SolverOptions& opts) :
    probes{ pIn },
    fes_{ fes }
{
    for (const auto& p: probes.exporterProbes) {
        exporterProbesCollection_.emplace(&p, buildParaviewDataCollectionInfo(p, fields));
    }

    for (const auto& p : probes.pointProbes) {
        pointProbesCollection_.emplace(&p, buildPointProbeCollectionInfo(p, fields));
    }

    for (const auto& p : probes.fieldProbes) {
        fieldProbesCollection_.emplace(&p, buildFieldProbeCollectionInfo(p, fields));
    }

    for (const auto& p : probes.nearFieldProbes) {
        auto dgfec{ dynamic_cast<const DG_FECollection*>(fes_.FEColl()) };
        if (!dgfec)
        {
            throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
        }
        auto reqs = std::make_unique<NearFieldReqs>(p, dgfec, fes_, fields);
        const bool local = reqs->hasLocalSurface();
        nearFieldReqs_.emplace(&p, std::move(reqs));
        {
            const MPI_Comm comm = getFESComm(fes_);
            const std::string parent_path =
                getSimulationCaseExportPath(caseName_) + "/NearToFarFieldProbes/" + p.name;
            if (isNodeRoot(comm)) {
                std::filesystem::create_directories(parent_path);
            }
            MPI_Barrier(comm);
        }
        if (local) {
            nearFieldProbesCollection_.emplace(&p, buildNearFieldDataCollectionInfo(p, fields));
        }
    }

    for (const auto& p: probes.domainSnapshotProbes) {
        domainSnapshotProbesCollection_.emplace(&p, buildDomainSnapshotDataCollection(p, fields));
    }
    
    finalTime_ = opts.final_time;
    fields_ = &fields;
    is_sgbc_solver_ = opts.is_sgbc_solver;
    preserve_existing_outputs_ = opts.resume_from_checkpoint;
}

static void fsyncWrittenFile(const std::filesystem::path& path)
{
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        return;
    }
    ::fsync(fd);
    ::close(fd);
}

void ProbesManager::flushOpenFiles()
{
    const auto root = std::filesystem::path(getSimulationCaseExportPath(caseName_));
    for (auto& entry : pointProbeFiles_) {
        if (!entry.second.is_open()) {
            continue;
        }
        entry.second.flush();
        fsyncWrittenFile(root / "PointProbes" / ("PointProbe" + std::to_string(entry.first) + ".dat"));
    }
    for (auto& entry : fieldProbeFiles_) {
        if (!entry.second.is_open()) {
            continue;
        }
        entry.second.flush();
        fsyncWrittenFile(root / "FieldProbes" / ("FieldProbe" + std::to_string(entry.first) + ".dat"));
    }
}

void ProbesManager::printTimingSummaryAndReset() const
{
#ifdef SHOW_TIMER_INFORMATION
    if (is_sgbc_solver_ || timingStats_.update_calls == 0) {
        return;
    }

    const int P = Mpi::WorldSize();
    double local_avg[7] = {
        timingStats_.exporter_ms / timingStats_.update_calls,
        timingStats_.field_ms    / timingStats_.update_calls,
        timingStats_.point_ms    / timingStats_.update_calls,
        timingStats_.nearfield_ms/ timingStats_.update_calls,
        timingStats_.snapshot_ms / timingStats_.update_calls,
        timingStats_.rcs_ms      / timingStats_.update_calls,
        timingStats_.mor_ms      / timingStats_.update_calls
    };

    std::vector<double> all_avg(7 * P);
    if (P > 1) {
        MPI_Gather(local_avg, 7, MPI_DOUBLE,
                   all_avg.data(), 7, MPI_DOUBLE,
                   0, getFESComm(fes_));
    } else {
        std::copy(local_avg, local_avg + 7, all_avg.data());
    }

    if (Mpi::WorldRank() == 0) {
        const bool defined[7] = {
            !probes.exporterProbes.empty(),
            !probes.fieldProbes.empty(),
            !probes.pointProbes.empty(),
            !probes.nearFieldProbes.empty(),
            !probes.domainSnapshotProbes.empty(),
            !probes.rcsSurfaceProbes.empty(),
            !probes.morStateProbes.empty()
        };
        const int widths[7] = {9, 9, 9, 10, 9, 9, 9};

        std::cout << "[Probe timing] avg of " << timingStats_.update_calls
                  << " updates, ms/update\n";
        std::cout << "  Rank  |  exporter     field     point  nearfield  snapshot       rcs       mor\n";
        std::cout << "  ------+-----------------------------------------------------------------------\n";
        for (int r = 0; r < P; ++r) {
            char cells[7][16];
            for (int c = 0; c < 7; ++c) {
                formatProbeTimingCell(
                    cells[c], sizeof(cells[c]), all_avg[r * 7 + c], defined[c], widths[c]);
            }
            char buf[256];
            std::snprintf(
                buf, sizeof(buf),
                "  %4d  | %s %s %s %s %s %s %s",
                r, cells[0], cells[1], cells[2], cells[3], cells[4], cells[5], cells[6]);
            std::cout << buf << '\n';
        }
        std::cout << std::flush;
    }

    timingStats_ = TimingStats{};
#endif
}

const FieldProbe& ProbesManager::getFieldProbe(const std::size_t i) const
{
    assert(i < probes.fieldProbes.size());
    return probes.fieldProbes[i];
}

const PointProbe& ProbesManager::getPointProbe(const std::size_t i) const
{
    assert(i < probes.pointProbes.size());
    return probes.pointProbes[i];
}

const ParGridFunction& getFieldView(const FieldProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields)
{
    switch (p.getFieldType()) {
    case FieldType::E:
        return fields.get(E, p.getDirection());
    case FieldType::H:
        return fields.get(H, p.getDirection());
    default:
        throw std::runtime_error("Invalid field type.");
    }
}

DenseMatrix pointVectorToDenseMatrixColumnVector(const Point& p)
{
    DenseMatrix r{(int)p.size(), 1 };
    for (auto i{ 0 }; i < p.size(); ++i) {
        r(i, 0) = p[i];
    }
    return r;
}

ProbesManager::PointProbeCollection
ProbesManager::buildPointProbeCollectionInfo(const PointProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{
    Array<int> elemIdArray;
    Array<IntegrationPoint> integPointArray;
    auto pointMatrix{ pointVectorToDenseMatrixColumnVector(p.getPoint()) };
    fes_.GetParMesh()->FindPoints(pointMatrix, elemIdArray, integPointArray);
    assert(elemIdArray.Size() == 1);
    assert(integPointArray.Size() == 1);
    FESPoint fesPoints { elemIdArray[0], integPointArray[0] };

    return { 
        fesPoints, 
        fields.get(E, X),
        fields.get(E, Y),
        fields.get(E, Z),
        fields.get(H, X),
        fields.get(H, Y),
        fields.get(H, Z)
    };
}

void ProbesManager::initPointFieldProbeExport()
{
    const MPI_Comm comm = getFESComm(fes_);

    if (probes.pointProbes.size()){
        auto base_path(getSimulationCaseExportPath(caseName_) + "/PointProbes/");
        
        if (cycle_ == 0 && !preserve_existing_outputs_) {
            if (isNodeRoot(comm)) {
                if (std::filesystem::exists(base_path)) {
                    std::filesystem::remove_all(base_path);
                }
                std::filesystem::create_directories(base_path);
            }
        }
        MPI_Barrier(comm);

        for (const auto& p : probes.pointProbes) {
            if(p.write){
                const auto& it{ pointProbesCollection_.find(&p) };
                if (it != pointProbesCollection_.end() && it->second.fesPoint.elementId != -2) {
                    std::ofstream myfile;
                    std::string path(base_path + "PointProbe" + std::to_string(p.getProbeID()) + ".dat");
                    std::vector<double> position = std::vector<double>({0.0, 0.0, 0.0});
                    for (auto i = 0; i < p.getPoint().size(); i++){
                        position[i] = p.getPoint()[i];
                    }
                    if (preserve_existing_outputs_) {
                        continue;
                    }
                    myfile.open(path, std::ios::trunc); 
                    if (myfile.is_open()) {
                        myfile << "PointProbe ID " << std::to_string(p.getProbeID()) << "\n";
                        myfile << "Spatial Position (X, Y, Z) \n";
                        myfile << std::scientific << std::setprecision(5);
                        myfile << std::to_string(position[0]) + " " + std::to_string(position[1]) + " " + std::to_string(position[2]) << "\n";
                        myfile << "Time (s) // Ex // Ey // Ez // Hx // Hy // Hz \n";
                    }
                    myfile.close();
                }
            }
        }
    }

    if (probes.fieldProbes.size()){
        auto base_path = (getSimulationCaseExportPath(caseName_) + "/FieldProbes/");
        
        if (cycle_ == 0 && !preserve_existing_outputs_) {
            if (isNodeRoot(comm)) {
                if (std::filesystem::exists(base_path)) {
                    std::filesystem::remove_all(base_path);
                }
                std::filesystem::create_directories(base_path);
            }
        }
        MPI_Barrier(comm); 

        for (const auto& p : probes.fieldProbes) {
            if(p.write){
                const auto& it{ fieldProbesCollection_.find(&p) };
                if (it != fieldProbesCollection_.end() && it->second.fesPoint.elementId != -2) {
                    std::ofstream myfile;
                    std::string path(base_path + "FieldProbe" + std::to_string(p.getProbeID()) + ".dat");
                    std::vector<double> position = std::vector<double>({0.0, 0.0, 0.0});
                    auto fieldpol = getFieldPolString(p.getFieldType(), p.getDirection());
                    for (auto i = 0; i < p.getPoint().size(); i++){
                        position[i] = p.getPoint()[i];
                    }
                    if (preserve_existing_outputs_) {
                        continue;
                    }
                    myfile.open(path, std::ios::trunc);
                    if (myfile.is_open()) {
                        myfile << "FieldProbe ID " << std::to_string(p.getProbeID()) << "\n";
                        myfile << "Spatial Position (X, Y, Z) \n";
                        myfile << std::scientific << std::setprecision(5);
                        myfile << std::to_string(position[0]) + " " + std::to_string(position[1]) + " " + std::to_string(position[2]) << "\n";
                        myfile << "Time (s) // " + fieldpol + "\n";
                    }
                    myfile.close();
                }
            }
        }
    }
}

ProbesManager::FieldProbeCollection
ProbesManager::buildFieldProbeCollectionInfo(const FieldProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{
    Array<int> elemIdArray;
    Array<IntegrationPoint> integPointArray;
    auto pointMatrix{ pointVectorToDenseMatrixColumnVector(p.getPoint()) };
    fes_.GetParMesh()->FindPoints(pointMatrix, elemIdArray, integPointArray);
    assert(elemIdArray.Size() == 1);
    assert(integPointArray.Size() == 1);
    FESPoint fesPoints{ elemIdArray[0], integPointArray[0] };

    return { fesPoints, getFieldView(p, fields) };
}

void isDGCollection(const FiniteElementSpace& fes)
{
    if (!dynamic_cast<const DG_FECollection*>(fes.FEColl()))
    {
        throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
    }
}

DataCollection ProbesManager::buildNearFieldDataCollectionInfo(
    const NearFieldProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& gFields) const
{
    isDGCollection(fes_);
    auto* reqs = nearFieldReqs_.at(&p).get();
    DataCollection res{ p.name, reqs->getSubMesh() };

    std::string path = getSimulationCaseExportPath(caseName_) + "/NearToFarFieldProbes/" + p.name
        + "/rank" + std::to_string(Mpi::WorldRank());
    std::filesystem::create_directories(path);
    res.SetPrefixPath(path);
    
    res.RegisterField("Ex.gf", &nearFieldReqs_.at(&p)->getConstField(E, X));
    res.RegisterField("Ey.gf", &nearFieldReqs_.at(&p)->getConstField(E, Y));
    res.RegisterField("Ez.gf", &nearFieldReqs_.at(&p)->getConstField(E, Z));
    res.RegisterField("Hx.gf", &nearFieldReqs_.at(&p)->getConstField(H, X));
    res.RegisterField("Hy.gf", &nearFieldReqs_.at(&p)->getConstField(H, Y));
    res.RegisterField("Hz.gf", &nearFieldReqs_.at(&p)->getConstField(H, Z));

    return res;
}

DomainSnapshotDataCollection ProbesManager::buildDomainSnapshotDataCollection(const DomainSnapshotProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const
{
    isDGCollection(fes_);
    DomainSnapshotDataCollection res(fes_, fields);
    return res;
}

void ProbesManager::setFinalTime(double final_time)
{
    finalTime_ = final_time;
    exporterContexts_.clear();
}

void ProbesManager::updateProbe(ExporterProbe& p, Time time)
{
    const MPI_Comm comm = getFESComm(fes_);
    int export_cycle = cycle_;

    if (p.save_every > 0.0 && finalTime_ > 0.0) {
        auto& ctx = exporterContexts_[&p];

        if (!ctx.initialized) {
            ctx.dt_save = p.save_every;
            ctx.next_save_time = 0.0;
            ctx.save_count = 0;
            ctx.finished = false;
            ctx.initialized = true;
        }

        if (ctx.finished) {
            return;
        }

        const double tol = ctx.dt_save * 1e-6;
        const double end_tol = std::max(tol, 1e-12);
        if (time < ctx.next_save_time - tol) {
            return;
        }

        export_cycle = ctx.save_count;
        ++ctx.save_count;

        if (ctx.next_save_time >= finalTime_ - end_tol) {
            ctx.finished = true;
        } else {
            const double next = ctx.next_save_time + ctx.dt_save;
            ctx.next_save_time = (next >= finalTime_ - end_tol) ? finalTime_ : next;
        }
    } else if (std::abs(time - finalTime_) >= 1e-8) {
        if (cycle_ % p.visSteps != 0) {
            return;
        }
    }

    auto it{ exporterProbesCollection_.find(&p) };
    assert(it != exporterProbesCollection_.end());
    auto& pd{ it->second };

    pd.SetCycle(export_cycle);
    pd.SetTime(std::round(time * 1e3) / 1e3);

    std::string base_dir = pd.GetPrefixPath() + pd.GetCollectionName();
    std::string cycle_dir = base_dir + "/Cycle" + to_padded_string(export_cycle, 6);
    
    if (isNodeRoot(comm)) {
        std::filesystem::create_directories(cycle_dir);
    }
    MPI_Barrier(comm);

#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (mfem::Device::Allows(mfem::Backend::CUDA) && fields_) {
        MFEM_STREAM_SYNC;
        fields_->allDOFs().HostRead();
        for (int d = X; d <= Z; ++d) {
            fields_->get(E, d).HostRead();
            fields_->get(H, d).HostRead();
        }
        fields_->get(E).HostRead();
        fields_->get(H).HostRead();
    }
#endif

    // ParaViewDataCollection::Save is per-rank I/O, not an MPI collective.
    // Without a post-Save barrier, a fast rank can leave step() and block in the
    // next stability/Mult Allreduce while a slow rank is still writing — hang.
    pd.Save();
    MPI_Barrier(comm);
}

void ProbesManager::updateProbe(FieldProbe& p, Time time)
{
    if (std::abs(time - finalTime_) >= 1e-8) {
        if (cycle_ % p.getVisSteps() != 0) {
            return;
        }
    }
    
    const auto& it{ fieldProbesCollection_.find(&p) };
    assert(it != fieldProbesCollection_.end());
    const auto& pC{ it->second };
    if (pC.fesPoint.elementId != -2){
#ifdef SEMBA_DGTD_ENABLE_CUDA
        if (mfem::Device::Allows(mfem::Backend::CUDA)) {
            const_cast<mfem::GridFunction&>(pC.field).HostRead();
        }
#endif
        real_t gf_value = pC.field.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);

        p.addFieldToMovies(time, gf_value);

        if(p.write){
            auto& myfile = fieldProbeFiles_[p.getProbeID()];
            if (!myfile.is_open()) {
                std::string path(getSimulationCaseExportPath(caseName_) + "/FieldProbes/" + "FieldProbe" + std::to_string(p.getProbeID()) + ".dat");
                myfile.open(path, std::ios::app);
            }
            if (myfile.is_open()) {
                myfile << std::scientific << std::setprecision(5);
                myfile << time / physicalConstants::speedOfLight_SI << " " << gf_value << "\n";
            }
        }
    }
}

void ProbesManager::updateProbe(PointProbe& p, Time time)
{
    if (std::abs(time - finalTime_) >= 1e-8) {
        if (cycle_ % p.getVisSteps() != 0) {
            return;
        }
    }
    
    const auto& it{ pointProbesCollection_.find(&p) };
    assert(it != pointProbesCollection_.end());
    const auto& pC{ it->second };
    if (pC.fesPoint.elementId != -2){
#ifdef SEMBA_DGTD_ENABLE_CUDA
        if (mfem::Device::Allows(mfem::Backend::CUDA)) {
            const_cast<mfem::GridFunction&>(pC.field_Ex).HostRead();
            const_cast<mfem::GridFunction&>(pC.field_Ey).HostRead();
            const_cast<mfem::GridFunction&>(pC.field_Ez).HostRead();
            const_cast<mfem::GridFunction&>(pC.field_Hx).HostRead();
            const_cast<mfem::GridFunction&>(pC.field_Hy).HostRead();
            const_cast<mfem::GridFunction&>(pC.field_Hz).HostRead();
        }
#endif
        FieldsForMovie f4FP;
        {
            f4FP.Ex = pC.field_Ex.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
            f4FP.Ey = pC.field_Ey.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
            f4FP.Ez = pC.field_Ez.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
            f4FP.Hx = pC.field_Hx.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
            f4FP.Hy = pC.field_Hy.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
            f4FP.Hz = pC.field_Hz.GetValue(pC.fesPoint.elementId, pC.fesPoint.iP);
        }
        p.addFieldsToMovies(time, f4FP);
        if(p.write){
            auto& myfile = pointProbeFiles_[p.getProbeID()];
            if (!myfile.is_open()) {
                std::string path(getSimulationCaseExportPath(caseName_) + "/PointProbes/" + "PointProbe" + std::to_string(p.getProbeID()) + ".dat");
                myfile.open(path, std::ios::app);
            }
            if (myfile.is_open()) {
                myfile << std::scientific << std::setprecision(5);
                myfile << time / physicalConstants::speedOfLight_SI << 
                " " << f4FP.Ex << " " << f4FP.Ey << " " << f4FP.Ez <<
                " " << f4FP.Hx << " " << f4FP.Hy << " " << f4FP.Hz << "\n";
            }
        }
    }
}

Fields<ParFiniteElementSpace, ParGridFunction> buildFieldsForProbe(const Fields<ParFiniteElementSpace, ParGridFunction>& src, ParFiniteElementSpace& fes)
{
    Fields<ParFiniteElementSpace, ParGridFunction> res(fes);
    for (auto f : { E, H }) {
        for (auto& d : { X, Y, Z }) {
            TransferMap tm(src.get(f, d), res.get(f, d));
            tm.Transfer(src.get(f, d), res.get(f, d));
        }
    }
    return res;
}

void ProbesManager::updateProbe(NearFieldProbe& p, Time time)
{
    if (std::abs(time - finalTime_) >= 1e-8) {
        if (cycle_ % p.expSteps != 0) {
            return;
        }
    }

    auto it{ nearFieldProbesCollection_.find(&p) };
    if (it == nearFieldProbesCollection_.end()) {
        return;
    }
    auto& dc{ it->second };
    dc.SetPrefixPath(getSimulationCaseExportPath(caseName_) + "/NearToFarFieldProbes/" + p.name + "/rank" + std::to_string(Mpi::WorldRank()));

    nearFieldReqs_.at(&p)->updateFields();

    dc.SetCycle(cycle_);
    dc.SetTime(time);
    dc.Save(); 

    std::string mesh_path{ dc.GetPrefixPath() + dc.GetCollectionName() + "/mesh" }; 
    auto mesh{ dc.GetMesh() };
    auto elemOrder = fes_.GetMesh()->GetElementTransformation(0)->Order();
    mesh->SetCurvature(elemOrder);
    mesh->Save(dc.GetPrefixPath() + "/mesh");

    std::string dir_name = dc.GetPrefixPath() + dc.GetCollectionName() + "_" + to_padded_string(dc.GetCycle(), 6) + "/time.txt";
    std::ofstream file;
    file.open(dir_name);
    file << time;
    file.close();
}

void ProbesManager::updateProbe(DomainSnapshotProbe& p, Time time)
{
    if (std::abs(time - finalTime_) >= 1e-8) {
        if (cycle_ % p.expSteps != 0) {
            return;
        }
    }

    auto it{ domainSnapshotProbesCollection_.find(&p) };
    assert(it != domainSnapshotProbesCollection_.end());
    auto& dc{ it->second };

#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (mfem::Device::Allows(mfem::Backend::CUDA) && fields_) {
        fields_->allDOFs().HostRead();
        for (int d = X; d <= Z; ++d) {
            fields_->get(E, d).HostRead();
            fields_->get(H, d).HostRead();
        }
    }
#endif

    std::string case_path = std::string(getSimulationCaseExportPath(caseName_) + "/DomainSnapshotProbes/");
    const MPI_Comm comm = getFESComm(fes_);
    
    if (cycle_ == 0) {
        if (isNodeRoot(comm)) {
            if (std::filesystem::exists(case_path)) {
                std::filesystem::remove_all(case_path);
            }
            std::filesystem::create_directories(case_path);
            std::filesystem::create_directories(case_path + "/meshes/");
        }
        MPI_Barrier(comm);

        dc.mesh.Save(case_path + "/meshes/mesh_rank" + std::to_string(Mpi::WorldRank()) , 16);
    }

    std::string rank_path = case_path + "/rank_" + std::to_string(Mpi::WorldRank());
    if (cycle_ == 0) {
        std::filesystem::create_directories(rank_path);
    }

    std::string folder_path = rank_path + "/cycle_" + to_padded_string(cycle_, 6) + "/";
    std::filesystem::create_directories(folder_path);

    dc.Save(folder_path);

    std::ofstream file(folder_path + "time.txt");
    file << time;
}

bool ProbesManager::needsHostSyncThisStep(Time time) const
{
    const bool at_end = std::abs(time - finalTime_) < 1e-8;

    auto dueBySteps = [&](int steps) {
        return at_end || (steps > 0 && (cycle_ % steps) == 0);
    };

    for (const auto& p : probes.exporterProbes) {
        if (p.save_every > 0.0 && finalTime_ > 0.0) {
            auto it = exporterContexts_.find(&p);
            if (it == exporterContexts_.end()) {
                return true; // first call initializes and may save t=0
            }
            const auto& ctx = it->second;
            if (ctx.finished) {
                continue;
            }
            const double tol = ctx.dt_save * 1e-6;
            if (time >= ctx.next_save_time - tol) {
                return true;
            }
        } else if (dueBySteps(p.visSteps)) {
            return true;
        }
    }

    for (const auto& p : probes.fieldProbes) {
        if (dueBySteps(p.getVisSteps())) {
            return true;
        }
    }
    for (const auto& p : probes.pointProbes) {
        if (dueBySteps(p.getVisSteps())) {
            return true;
        }
    }
    for (const auto& p : probes.nearFieldProbes) {
        if (dueBySteps(p.expSteps)) {
            return true;
        }
    }
    for (const auto& p : probes.domainSnapshotProbes) {
        if (dueBySteps(p.expSteps)) {
            return true;
        }
    }
    for (const auto& p : probes.rcsSurfaceProbes) {
        if (dueBySteps(p.expSteps)) {
            return true;
        }
    }
    for (const auto& p : probes.morStateProbes) {
        if (p.saves <= 0) {
            continue;
        }
        if (time < p.record_time_start - 1e-12 || time > p.record_time_final + 1e-12) {
            continue;
        }
        auto it = morStateContexts_.find(&p);
        if (it == morStateContexts_.end()) {
            return true;
        }
        const auto& ctx = it->second;
        if (ctx.save_count >= p.saves) {
            continue;
        }
        const double tol = (ctx.dt_save > 0.0) ? ctx.dt_save * 1e-6 : 1e-12;
        if (time >= ctx.next_save_time - tol) {
            return true;
        }
    }

    return false;
}

void ProbesManager::updateProbes(Time t)
{
#ifdef SHOW_TIMER_INFORMATION
    using clock = std::chrono::steady_clock;
    auto start = clock::now();
#endif
    for (auto& p : probes.exporterProbes) {
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_exporter = clock::now();
#endif
    for (auto& p : probes.fieldProbes) {
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_field = clock::now();
#endif
    for (auto& p : probes.pointProbes) {
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_point = clock::now();
#endif
    for (auto& p : probes.nearFieldProbes) {
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_nearfield = clock::now();
#endif
    for (auto& p : probes.domainSnapshotProbes){
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_snapshot = clock::now();
#endif
    for (auto& p : probes.rcsSurfaceProbes) {
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_rcs = clock::now();
#endif
    for (auto& p : probes.morStateProbes) {
        updateProbe(p, t);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_mor = clock::now();
    timingStats_.exporter_ms += std::chrono::duration<double, std::milli>(after_exporter - start).count();
    timingStats_.field_ms    += std::chrono::duration<double, std::milli>(after_field - after_exporter).count();
    timingStats_.point_ms    += std::chrono::duration<double, std::milli>(after_point - after_field).count();
    timingStats_.nearfield_ms+= std::chrono::duration<double, std::milli>(after_nearfield - after_point).count();
    timingStats_.snapshot_ms += std::chrono::duration<double, std::milli>(after_snapshot - after_nearfield).count();
    timingStats_.rcs_ms      += std::chrono::duration<double, std::milli>(after_rcs - after_snapshot).count();
    timingStats_.mor_ms      += std::chrono::duration<double, std::milli>(after_mor - after_rcs).count();
    timingStats_.update_calls++;
#endif

    cycle_++;
}

void NearFieldReqs::updateFields()
{
    if (!tMaps_) {
        return;
    }
    tMaps_->transferFields(gFields_, *fields_);
}

NearFieldReqs::NearFieldReqs(
    const NearFieldProbe& p, const DG_FECollection* fec, ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>& global) :
    ntff_smsh_{ NearToFarFieldSubMesher(*fes.GetMesh(), fes, buildSurfaceMarker(p.tags, fes)) },
    gFields_{ global }
{
    if (!ntff_smsh_.hasLocalSurface()) {
        return;
    }
    sfes_ = std::make_unique<FiniteElementSpace>(ntff_smsh_.getSubMesh(), fec);
    fields_ = std::make_unique<Fields<FiniteElementSpace, GridFunction>>(*sfes_);
    tMaps_ = std::make_unique<TransferMaps>(gFields_, *fields_);
    updateFields();
}

void ProbesManager::updateProbe(RCSSurfaceProbe& p, Time t)
{
    auto it = rcsSurfaceExporters_.find(&p);
    if (it != rcsSurfaceExporters_.end()) {
        it->second->write(t, cycle_, finalTime_);
    }
}

void ProbesManager::initRCSSurfaceExporters()
{
    if (!fields_) return;
    for (const auto& p : probes.rcsSurfaceProbes) {
        auto dgfec = dynamic_cast<const DG_FECollection*>(fes_.FEColl());
        if (!dgfec) {
            throw std::runtime_error("The FiniteElementCollection in the FiniteElementSpace is not DG.");
        }
        rcsSurfaceExporters_.emplace(&p,
            std::make_unique<RCSSurfaceExporter>(
                p, dgfec, fes_, *fields_, caseName_, preserve_existing_outputs_));
    }
}

void ProbesManager::recalculateExportSteps(double dt)
{
    if (dt <= 0.0) return;

    const int totalSteps = static_cast<int>(std::ceil(finalTime_ / dt));

    auto stepsFromSaves = [&](int saves) -> int {
        if (saves <= 0) return 1;
        return std::max(1, totalSteps / saves);
    };

    // ExporterProbe::save_every is absolute time; no step remapping needed.
    for (auto& p : probes.nearFieldProbes) {
        if (p.saves > 0) p.expSteps = stepsFromSaves(p.saves);
    }
    for (auto& p : probes.rcsSurfaceProbes) {
        if (p.saves > 0) p.expSteps = stepsFromSaves(p.saves);
    }
    for (auto& p : probes.domainSnapshotProbes) {
        if (p.saves > 0) p.expSteps = stepsFromSaves(p.saves);
    }
    for (auto& p : probes.fieldProbes) {
        if (p.getSaves() > 0) p.setVisSteps(stepsFromSaves(p.getSaves()));
    }
    for (auto& p : probes.pointProbes) {
        if (p.getSaves() > 0) p.setVisSteps(stepsFromSaves(p.getSaves()));
    }
}

void ProbesManager::updateProbe(MORStateProbe& p, Time time)
{
    if (p.saves <= 0) return;
    if (time < p.record_time_start - 1e-12) return;
    if (time > p.record_time_final + 1e-12) return;

    auto& ctx = morStateContexts_[&p];

    if (!ctx.initialized) {
        const MPI_Comm comm = getFESComm(fes_);
        ctx.dt_save = (p.saves > 1)
            ? (p.record_time_final - p.record_time_start) / (p.saves - 1)
            : 0.0;
        ctx.next_save_time = p.record_time_start;
        ctx.save_count = 0;
        ctx.export_dir = getSimulationCaseExportPath(caseName_) + "/MORStateProbes/" + p.name;

        if (isNodeRoot(comm)) {
            std::filesystem::create_directories(ctx.export_dir);
        }
        MPI_Barrier(comm);

        ctx.initialized = true;
    }

    if (ctx.save_count >= p.saves) return;

    double tol = (ctx.dt_save > 0.0) ? ctx.dt_save * 1e-6 : 1e-12;
    if (time < ctx.next_save_time - tol) return;

    // Export E/H state only (first 6×ndofs block; ψ is internal).
    const auto& all_dofs = fields_->allDOFs();
#ifdef SEMBA_DGTD_ENABLE_CUDA
    if (mfem::Device::Allows(mfem::Backend::CUDA)) {
        const_cast<mfem::Vector&>(all_dofs).HostRead();
    }
#endif
    const int export_size = fields_->fieldBlockSize();
    std::string file_path = ctx.export_dir + "/x_" + std::to_string(ctx.save_count);

    std::ofstream ofs(file_path);
    if (ofs.is_open()) {
        ofs << std::scientific << std::setprecision(16);
        ofs << time << "\n";
        ofs << export_size << "\n";
        for (int i = 0; i < export_size; ++i) {
            ofs << all_dofs[i] << "\n";
        }
        ofs.close();
    }

    // Export TFSF source function u if available.
    // Use host-only eval so we never HostRead/Read the device-resident
    // cached_tfsf_fields_ (that poisons the next GPU Mult TFSF eval).
    if (srcmngr_ && tfsf_mapping_ && srcmngr_->hasDirectEval()) {
        std::vector<double> u_combined;
        srcmngr_->evalTimeVarFieldDirectToHost(time, u_combined);
        const int tfsf_total_size = static_cast<int>(u_combined.size());

        std::string u_path = ctx.export_dir + "/u_" + std::to_string(ctx.save_count);
        std::ofstream ofs_u(u_path);
        if (ofs_u.is_open()) {
            ofs_u << std::scientific << std::setprecision(16);
            ofs_u << time << "\n";
            ofs_u << tfsf_total_size << "\n";
            for (int i = 0; i < tfsf_total_size; ++i) {
                ofs_u << u_combined[static_cast<std::size_t>(i)] << "\n";
            }
            ofs_u.close();
        }
    }

    ctx.save_count++;
    if (p.saves > 1) {
        ctx.next_save_time = p.record_time_start + ctx.save_count * ctx.dt_save;
    } else {
        ctx.next_save_time = p.record_time_final + 1.0;
    }
}

int expectedProbeSamples(int saved_cycle, int vis_steps, bool at_final)
{
    if (saved_cycle <= 0) {
        return 0;
    }
    if (vis_steps <= 0) {
        return at_final ? 1 : 0;
    }
    int count = (saved_cycle - 1) / vis_steps + 1;
    if (at_final && ((saved_cycle - 1) % vis_steps) != 0) {
        ++count;
    }
    return count;
}

int cycleInPvdLine(const std::string& line)
{
    const auto pos = line.find("Cycle");
    if (pos == std::string::npos) {
        return -1;
    }
    std::size_t i = pos + 5;
    if (i >= line.size() || line[i] < '0' || line[i] > '9') {
        return -1;
    }
    int value = 0;
    while (i < line.size() && line[i] >= '0' && line[i] <= '9') {
        value = value * 10 + (line[i] - '0');
        ++i;
    }
    return value;
}

namespace {

void truncateDataLines(const std::filesystem::path& path, int header_lines, int keep_data_lines)
{
    if (!std::filesystem::exists(path)) {
        return;
    }
    std::ifstream in(path);
    if (!in) {
        return;
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    in.close();
    const int keep = std::min(
        header_lines + std::max(keep_data_lines, 0),
        static_cast<int>(lines.size()));
    if (keep >= static_cast<int>(lines.size())) {
        return;
    }
    std::ofstream out(path, std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Cannot trim " + path.string());
    }
    for (int i = 0; i < keep; ++i) {
        out << lines[i] << '\n';
    }
}

int trailingIndex(const std::string& name, const std::string& prefix)
{
    if (name.size() <= prefix.size() || name.compare(0, prefix.size(), prefix) != 0) {
        return -1;
    }
    int value = 0;
    for (std::size_t i = prefix.size(); i < name.size(); ++i) {
        if (name[i] < '0' || name[i] > '9') {
            return -1;
        }
        value = value * 10 + (name[i] - '0');
    }
    return value;
}

void removeIndexedDirs(const std::filesystem::path& parent, const std::string& prefix, int drop_from)
{
    if (!std::filesystem::exists(parent)) {
        return;
    }
    for (const auto& entry : std::filesystem::directory_iterator(parent)) {
        if (!entry.is_directory()) {
            continue;
        }
        const int index = trailingIndex(entry.path().filename().string(), prefix);
        if (index >= drop_from) {
            std::filesystem::remove_all(entry.path());
        }
    }
}

void removeIndexedFiles(const std::filesystem::path& dir, const std::string& prefix, int drop_from)
{
    if (!std::filesystem::exists(dir)) {
        return;
    }
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const int index = trailingIndex(entry.path().filename().string(), prefix);
        if (index >= drop_from) {
            std::filesystem::remove(entry.path());
        }
    }
}

void trimPvd(const std::filesystem::path& path, int drop_from)
{
    if (!std::filesystem::exists(path)) {
        return;
    }
    std::ifstream in(path);
    if (!in) {
        return;
    }
    std::vector<std::string> kept;
    std::string line;
    bool dropped = false;
    while (std::getline(in, line)) {
        const int cycle = cycleInPvdLine(line);
        if (cycle >= drop_from) {
            dropped = true;
            continue;
        }
        kept.push_back(line);
    }
    in.close();
    if (!dropped) {
        return;
    }
    std::ofstream out(path, std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Cannot trim " + path.string());
    }
    for (const auto& kept_line : kept) {
        out << kept_line << '\n';
    }
}

}  // namespace

std::vector<ExporterCursor> ProbesManager::captureExporterCursors() const
{
    std::vector<ExporterCursor> out;
    out.reserve(probes.exporterProbes.size());
    for (const auto& p : probes.exporterProbes) {
        ExporterCursor cursor;
        cursor.name = p.name;
        const auto it = exporterContexts_.find(&p);
        if (it != exporterContexts_.end()) {
            cursor.save_count = it->second.save_count;
            cursor.next_save_time = it->second.next_save_time;
            cursor.dt_save = it->second.dt_save;
            cursor.initialized = it->second.initialized;
            cursor.finished = it->second.finished;
        }
        out.push_back(std::move(cursor));
    }
    return out;
}

std::vector<MorCursor> ProbesManager::captureMorCursors() const
{
    std::vector<MorCursor> out;
    out.reserve(probes.morStateProbes.size());
    for (const auto& p : probes.morStateProbes) {
        MorCursor cursor;
        cursor.name = p.name;
        const auto it = morStateContexts_.find(&p);
        if (it != morStateContexts_.end()) {
            cursor.save_count = it->second.save_count;
            cursor.next_save_time = it->second.next_save_time;
            cursor.dt_save = it->second.dt_save;
            cursor.initialized = it->second.initialized;
        }
        out.push_back(std::move(cursor));
    }
    return out;
}

void ProbesManager::restoreCheckpointCursors(
    int cycle,
    const std::vector<ExporterCursor>& exporters,
    const std::vector<MorCursor>& mor)
{
    if (exporters.size() != probes.exporterProbes.size()) {
        throw std::runtime_error("Checkpoint exporter cursors do not match this JSON.");
    }
    if (mor.size() != probes.morStateProbes.size()) {
        throw std::runtime_error("Checkpoint MOR cursors do not match this JSON.");
    }
    cycle_ = cycle;
    for (std::size_t i = 0; i < exporters.size(); ++i) {
        if (exporters[i].name != probes.exporterProbes[i].name) {
            throw std::runtime_error("Checkpoint exporter name does not match this JSON.");
        }
        if (!exporters[i].initialized) {
            continue;
        }
        auto& ctx = exporterContexts_[&probes.exporterProbes[i]];
        ctx.save_count = exporters[i].save_count;
        ctx.next_save_time = exporters[i].next_save_time;
        ctx.dt_save = exporters[i].dt_save;
        ctx.initialized = true;
        ctx.finished = exporters[i].finished;
    }
    for (std::size_t i = 0; i < mor.size(); ++i) {
        if (mor[i].name != probes.morStateProbes[i].name) {
            throw std::runtime_error("Checkpoint MOR probe name does not match this JSON.");
        }
        if (!mor[i].initialized) {
            continue;
        }
        auto& ctx = morStateContexts_[&probes.morStateProbes[i]];
        ctx.save_count = mor[i].save_count;
        ctx.next_save_time = mor[i].next_save_time;
        ctx.dt_save = mor[i].dt_save;
        ctx.initialized = true;
        ctx.export_dir = getSimulationCaseExportPath(caseName_) + "/MORStateProbes/" + probes.morStateProbes[i].name;
    }
    for (auto& entry : exporterProbesCollection_) {
        entry.second.UseRestartMode(true);
    }
}

void ProbesManager::trimProbeOutput(double time)
{
    const bool at_final = std::abs(time - finalTime_) < 1e-8;
    const int saved_cycle = cycle_;

    for (auto& entry : pointProbeFiles_) {
        if (entry.second.is_open()) {
            entry.second.close();
        }
    }
    for (auto& entry : fieldProbeFiles_) {
        if (entry.second.is_open()) {
            entry.second.close();
        }
    }

    const auto case_root = std::filesystem::path(getSimulationCaseExportPath(caseName_));

    for (const auto& p : probes.pointProbes) {
        if (!p.write) {
            continue;
        }
        const auto it = pointProbesCollection_.find(&p);
        if (it == pointProbesCollection_.end() || it->second.fesPoint.elementId == -2) {
            continue;
        }
        truncateDataLines(
            case_root / "PointProbes" / ("PointProbe" + std::to_string(p.getProbeID()) + ".dat"),
            4,
            expectedProbeSamples(saved_cycle, p.getVisSteps(), at_final));
    }
    for (const auto& p : probes.fieldProbes) {
        if (!p.write) {
            continue;
        }
        const auto it = fieldProbesCollection_.find(&p);
        if (it == fieldProbesCollection_.end() || it->second.fesPoint.elementId == -2) {
            continue;
        }
        truncateDataLines(
            case_root / "FieldProbes" / ("FieldProbe" + std::to_string(p.getProbeID()) + ".dat"),
            4,
            expectedProbeSamples(saved_cycle, p.getVisSteps(), at_final));
    }

    const int rank = Mpi::WorldRank();
    for (const auto& p : probes.nearFieldProbes) {
        const auto rank_dir = case_root / "NearToFarFieldProbes" / p.name
            / ("rank" + std::to_string(rank));
        removeIndexedDirs(rank_dir, p.name + "_", saved_cycle);
    }
    for (const auto& p : probes.domainSnapshotProbes) {
        (void)p;
        const auto rank_dir = case_root / "DomainSnapshotProbes" / ("rank_" + std::to_string(rank));
        removeIndexedDirs(rank_dir, "cycle_", saved_cycle);
    }
    for (auto& entry : rcsSurfaceExporters_) {
        entry.second->truncateSnapshotsAfter(time);
    }

    if (rank == 0) {
        const auto paraview_root = std::filesystem::path("exports/ParaView") / getRunModeTag();
        for (const auto& p : probes.exporterProbes) {
            int drop_from = saved_cycle;
            if (p.save_every > 0.0) {
                const auto it = exporterContexts_.find(&p);
                drop_from = (it == exporterContexts_.end()) ? 0 : it->second.save_count;
            }
            const auto parent = paraview_root / p.name;
            removeIndexedDirs(parent, "Cycle", drop_from);
            trimPvd(parent / (p.name + ".pvd"), drop_from);
        }
        for (const auto& p : probes.morStateProbes) {
            const auto it = morStateContexts_.find(&p);
            const int drop_from = (it == morStateContexts_.end() || !it->second.initialized)
                ? 0 : it->second.save_count;
            const auto dir = case_root / "MORStateProbes" / p.name;
            removeIndexedFiles(dir, "x_", drop_from);
            removeIndexedFiles(dir, "u_", drop_from);
        }
    }
}

}