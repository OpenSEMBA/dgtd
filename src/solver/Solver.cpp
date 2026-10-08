#include "Solver.h"
#include "Checkpoint.h"
#include "components/SCPMLLayout.h"
#include "evolution/EvolutionOptions.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <fcntl.h>
#include <atomic>
#include <csignal>
#include <cstring>
#include <map>
#include <optional>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#ifdef SEMBA_DGTD_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace maxwell {

size_t getCurrentMemoryUsage();
size_t getPeakMemoryUsage();

void Solver::flushProbeFiles()
{
    probesManager_.flushOpenFiles();
}

void Solver::sampleInitializationMemory()
{
    const auto cur = getCurrentMemoryUsage();
    if (cur > initMemPeakSampled_) {
        initMemPeakSampled_ = cur;
    }
}

void Solver::sampleTemporalMemory()
{
    const auto cur = getCurrentMemoryUsage();
    temporalMemSumSampled_ += static_cast<long double>(cur);
    temporalMemSampleCount_++;
    if (cur > temporalMemPeakSampled_) {
        temporalMemPeakSampled_ = cur;
    }
}

void warnGaussianPulseVsMesh(
	mfem::ParMesh& mesh,
	int order,
	const Sources& sources,
	int comm_rank)
{
	double min_spread = std::numeric_limits<double>::infinity();
	double max_carrier_hz = 0.0;
	const double c_si = physicalConstants::speedOfLight_SI;

	for (const auto& src : sources) {
		if (const auto* gap = dynamic_cast<const DeltaGapSource*>(src.get())) {
			min_spread = std::min(min_spread, gap->spread());
			continue;
		}
		const auto* tf = dynamic_cast<const TotalField*>(src.get());
		if (tf == nullptr) {
			continue;
		}
		const EHFieldFunction* eh = tf->function();
		if (const auto* pw = dynamic_cast<const Planewave*>(eh)) {
			if (const auto* g = dynamic_cast<const Gaussian*>(pw->function())) {
				min_spread = std::min(min_spread, g->spread());
			}
			else if (const auto* mg = dynamic_cast<const ModulatedGaussian*>(pw->function())) {
				min_spread = std::min(min_spread, mg->spread());
				max_carrier_hz = std::max(max_carrier_hz, mg->frequency() * c_si);
			}
		}
		else if (const auto* coax = dynamic_cast<const CoaxialMode*>(eh)) {
			min_spread = std::min(min_spread, coax->spread());
		}
		else if (const auto* dip = dynamic_cast<const DerivGaussDipole*>(eh)) {
			min_spread = std::min(min_spread, dip->spread());
		}
	}

	if (!(min_spread > 0.0) || !std::isfinite(min_spread)) {
		return;
	}

	double local_sum = 0.0;
	double local_min = std::numeric_limits<double>::infinity();
	const int local_ne = mesh.GetNE();
	for (int e = 0; e < local_ne; ++e) {
		const double h = mesh.GetElementSize(e);
		local_sum += h;
		local_min = std::min(local_min, h);
	}

	int global_ne = 0;
	double global_sum = 0.0;
	double global_min = 0.0;
	MPI_Comm comm = mesh.GetComm();
	MPI_Allreduce(&local_ne, &global_ne, 1, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
	if (global_ne <= 0 || !(global_sum > 0.0) || !std::isfinite(global_min)) {
		return;
	}

	const double h_avg = global_sum / static_cast<double>(global_ne);
	const double n_ppw = std::max(3.0, 15.0 / static_cast<double>(order + 1));
	const double f_mesh_hz = c_si / (n_ppw * h_avg);
	const double f_1e_hz = c_si / (2.0 * M_PI * min_spread) + max_carrier_hz;
	if (!(f_1e_hz > f_mesh_hz)) {
		return;
	}

	const double spread_ok = c_si / (2.0 * M_PI * std::max(f_mesh_hz - max_carrier_hz, 1.0));
	if (comm_rank != 0) {
		return;
	}

	std::cerr << std::setprecision(4)
	          << "Warning: Gaussian 1/e power sits at " << f_1e_hz / 1e6
	          << " MHz (spread=" << min_spread << " light-metres)"
	          << (max_carrier_hz > 0.0 ? " including carrier" : "")
	          << ".\n  Mesh order p=" << order
	          << ", mean h=" << h_avg << " m, min h=" << global_min << " m.\n"
	          << "  Estimated resolved band ~ " << f_mesh_hz / 1e6
	          << " MHz (~" << n_ppw << " elements/λ at mean h, order "
	          << order << ").\n"
	          << "  Consider f_1e <= " << f_mesh_hz
	          << " Hz or spread >= " << spread_ok
	          << ", or refine the mesh. Continuing.\n";
}

Solver::~Solver() = default; 

std::unique_ptr<mfem::ParFiniteElementSpace> buildFiniteElementSpace(mfem::ParMesh* m, mfem::FiniteElementCollection* fec)
{
    auto fes = std::make_unique<mfem::ParFiniteElementSpace>(m, fec);
    fes->ExchangeFaceNbrData();
    return fes;
}

std::unique_ptr<mfem::TimeDependentOperator> Solver::assignEvolutionOperator()
{
    if (opts_.evolution.op == EvolutionOperatorType::Hesthaven) {
        return std::make_unique<HesthavenEvolution>(*fes_, model_, sourcesManager_, opts_.evolution);
    }
    else if (opts_.evolution.op == EvolutionOperatorType::Global) {
        auto global_evol = std::make_unique<GlobalEvolution>(*fes_, model_, sourcesManager_, opts_.evolution, probesManager_.probes, opts_.final_time);
        globalEvol_cache_ = global_evol.get();  // Cache the pointer
        return global_evol;
    }
    else if (opts_.evolution.op == EvolutionOperatorType::Maxwell){
        ProblemDescription pd(model_, probesManager_.probes, sourcesManager_.sources, opts_.evolution);
        return std::make_unique<MaxwellEvolution>(pd, *fes_, sourcesManager_);
    }
    else {
        throw std::runtime_error("Unknown evolution operator type.");
    }
}

void Solver::assignODESolver()
{
    switch (static_cast<ode_type>(opts_.ode_type))
    {
        case ode_type::RK4:
            odeSolver_ = std::make_unique<mfem::RK4Solver>();
            break;

        case ode_type::BackwardEuler:
            odeSolver_ = std::make_unique<mfem::BackwardEulerSolver>();
            break;

        case ode_type::Trapezoidal:
            odeSolver_ = std::make_unique<mfem::TrapezoidalRuleSolver>();
            break;

        case ode_type::ImplicitMidpoint:
            odeSolver_ = std::make_unique<mfem::ImplicitMidpointSolver>();
            break;

        case ode_type::SDIRK33:
            odeSolver_ = std::make_unique<mfem::SDIRK33Solver>();
            break;

        case ode_type::SDIRK23:  // L-stable. Parent PML/Debye/Lorentz/SGBC stays on RK4.
            odeSolver_ = std::make_unique<mfem::SDIRK23Solver>(/*gamma_opt=*/2);
            break;

        case ode_type::SDIRK34:
            odeSolver_ = std::make_unique<mfem::SDIRK34Solver>();
            break;
            
        default:
            throw std::runtime_error(
                "Wrong ode type defined in json. See ode_type in SolverOptions for available inputs.");
    }
}

Solver::Solver(
    const Model& model,
    const Probes& probes,
    const Sources& sources,
    const SolverOptions& options) :
    opts_{ options },
    model_{ model },
    fec_{ opts_.evolution.order, model_.getMesh().Dimension(), opts_.basis_type},
    fes_{ buildFiniteElementSpace(& model_.getMesh(), &fec_) },
    fields_{ *fes_,
             computePMLAuxSize(
                 model.getPMLProperties(),
                 fes_->GetNDofs(),
                 fes_->GetMesh()->Dimension())
             + model.debyeAuxSize(fes_->GetNDofs())
             + model.lorentzAuxSize(fes_->GetNDofs()) },
    sourcesManager_{ sources, *fes_, fields_ },
    probesManager_ { probes , *fes_, fields_, opts_ },
    time_{0.0}
{
    auto initStartTime = std::chrono::steady_clock::now();
    initMemBaseline_ = getCurrentMemoryUsage();
    initMemPeakSampled_ = initMemBaseline_;

	MPI_Comm comm = model_.getMesh().GetComm();
    int comm_rank;
    MPI_Comm_rank(comm, &comm_rank);

    if (comm_rank == 0){
        checkOptionsAreValid(opts_);
    }
    warnGaussianPulseVsMesh(model_.getMesh(), opts_.evolution.order, sources, comm_rank);

    const bool dispersive = model_.hasDebye() || model_.hasLorentz();
    if (dispersive && opts_.evolution.op != EvolutionOperatorType::Global) {
        throw std::runtime_error(
            "Debye and Lorentz materials require evolution_operator \"global\".");
    }
    const bool implicit_blocked =
        model_.hasPML() || dispersive || !model_.getSGBCProperties().empty();
    if (implicit_blocked && opts_.ode_type != ode_type::RK4) {
        throw std::runtime_error(
            "Implicit ODE integrators do not support Cartesian PML, Debye, Lorentz, or SGBC. Use ode_type RK4.");
    }
    if (dispersive && opts_.evolution.spectral) {
        throw std::runtime_error(
            "Spectral analysis does not include Debye or Lorentz polarization.");
    }

    if (opts_.evolution.spectral == true) {
        performSpectralAnalysis(*fes_.get(), model_, opts_.evolution);
    }
    sampleInitializationMemory();

    if (model_.hasPML()) {
        model_.initializePMLProfiles(comm_rank, opts_.evolution.order);
    }

    evolTDO_ = assignEvolutionOperator();
    evolTDO_->SetTime(time_);
    
    // Initialize TFSF export in ProbesManager if TFSF is present and MOR state probes exist
    if (globalEvol_cache_ && !probesManager_.probes.morStateProbes.empty()) {
        probesManager_.initTFSFExport(&sourcesManager_, &(globalEvol_cache_->getTFSFMapping()));
    }
    
    sampleInitializationMemory();

    if (opts_.time_step == 0.0) {
        dt_ = estimateTimeStep(model_, opts_, *fes_, evolTDO_.get());
    }
    else {
        dt_ = opts_.time_step;
    }

    // Now that the real time step is known, ensure all probe export intervals
    // that were specified via "saves" are consistent with the actual dt.
    // This is especially important when using the automatic time stepper
    // (time_step == 0.0), where the interval could not be calculated correctly
    // at parse time.
    probesManager_.recalculateExportSteps(dt_);

    assignODESolver();
    odeSolver_->Init(*evolTDO_);
    sampleInitializationMemory();

    probesManager_.setCaseName(model_.meshName_);
    probesManager_.initPointFieldProbeExport();
    if (!opts_.resume_from_checkpoint) {
        probesManager_.updateProbes(time_);
    }
    sampleInitializationMemory();

    auto initEndTime = std::chrono::steady_clock::now();

    if (!opts_.is_sgbc_solver && !opts_.resume_from_checkpoint) {
        std::filesystem::path simExpPath(getSimulationCaseExportPath(model.meshName_) + "/SimulationStats/");
        std::string path = (simExpPath / ("statistics_rank" + std::to_string(comm_rank) + ".dat")).string();

        std::ofstream myfile(path, std::ios::app);
        if (myfile.is_open()) {
            auto runtime = std::chrono::duration<double>(initEndTime - initStartTime).count();
            myfile << std::scientific << std::setprecision(5);
            myfile << "Initialization Time: " << runtime << " (s)\n";
            myfile << std::defaultfloat;
            myfile << "Initialization Baseline Memory (B): " << initMemBaseline_ << "\n";
            myfile << "Initialization Peak Sampled Memory (B): " << initMemPeakSampled_ << "\n";
            myfile << "Initialization Peak Increase (B): " << (initMemPeakSampled_ - initMemBaseline_) << "\n";
            myfile << "Assembly Peak Memory Consumption (B): " << getPeakMemoryUsage() << "\n";
            myfile.close();
        } else {
            std::cerr << "Rank " << comm_rank << " failed to open file: " << path << "\n";
        }
    }

    if (!opts_.checkpoint_json_path.empty()) {
        enableCheckpoint(opts_.checkpoint_json_path, opts_.checkpoint_mesh_path);
    }
    if (!opts_.checkpoint_directory.empty()) {
        restoreCheckpoint(opts_.checkpoint_directory);
    }
}

void Solver::checkOptionsAreValid(const SolverOptions& opts) const
{
    
    if ((opts.evolution.order < 0) ||
        (opts.final_time < 0)) {
        throw std::runtime_error("Incorrect parameters in Options");
    }

    if (!(opts.checkpoint_percent >= 0.0 && opts.checkpoint_percent <= 100.0)) {
        throw std::runtime_error(
            "solver_options.checkpoint_percent must be between 0 and 100.");
    }

    if (opts.cfl <= 0.0) {
        throw std::runtime_error("CFL must be positive");
    }

}

const PointProbe& Solver::getPointProbe(const std::size_t probe) const 
{ 
    return probesManager_.getPointProbe(probe); 
}

const FieldProbe& Solver::getFieldProbe(const std::size_t probe) const
{
    return probesManager_.getFieldProbe(probe);
}

double getMinimumInterNodeDistance(FiniteElementSpace& fes)
{
    GridFunction nodes(&fes);
    fes.GetMesh()->GetNodes(nodes);
    double res{ std::numeric_limits<double>::max() };
    for (int e = 0; e < fes.GetMesh()->ElementToElementTable().Size(); ++e) {
        Array<int> dofs;
        fes.GetElementDofs(e, dofs);
        if (dofs.Size() == 1) {
            res = std::min(res, fes.GetMesh()->GetElementSize(e));
        }
        else {
            for (int i = 0; i < dofs.Size(); ++i) {
                for (int j = i + 1; j < dofs.Size(); ++j) {
                    res = std::min(res, std::abs(nodes[dofs[i]] - nodes[dofs[j]]));
                }
            }
        }
    }
    return res;
}

bool checkIfElemTypeInMesh(const mfem::Mesh& mesh, const mfem::Element::Type& type)
{
    for (int e = 0; e < mesh.GetNE(); ++e) {
        if (mesh.GetElementType(e) == type) {
            return true;
        }
    }
    return false;
}

mfem::Vector getTimeStepScale(mfem::Mesh& mesh)
{
    const int dim = mesh.Dimension();
    const int sdim = mesh.SpaceDimension();
    mfem::Vector vol(mesh.GetNE()), dtscale(mesh.GetNE());
    for (int e = 0; e < mesh.GetNE(); ++e) {
        double faceAreaSum = 0.0;

        mfem::Array<int> elem_faces, ori;
        if (dim == 2) {
            mesh.GetElementEdges(e, elem_faces, ori);
        } else {
            mesh.GetElementFaces(e, elem_faces, ori);
        }

        for (int fi = 0; fi < elem_faces.Size(); ++fi) {
            if (dim == 2) {
                // Edge length from vertices
                mfem::Array<int> v;
                mesh.GetEdgeVertices(elem_faces[fi], v);
                double len = 0.0;
                for (int d = 0; d < sdim; ++d) {
                    double diff = mesh.GetVertex(v[1])[d] - mesh.GetVertex(v[0])[d];
                    len += diff * diff;
                }
                faceAreaSum += std::sqrt(len);
            } else {
                // Triangle area from vertices
                mfem::Array<int> v;
                mesh.GetFaceVertices(elem_faces[fi], v);
                double e1[3], e2[3];
                for (int d = 0; d < 3; ++d) {
                    e1[d] = mesh.GetVertex(v[1])[d] - mesh.GetVertex(v[0])[d];
                    e2[d] = mesh.GetVertex(v[2])[d] - mesh.GetVertex(v[0])[d];
                }
                double cx = e1[1]*e2[2] - e1[2]*e2[1];
                double cy = e1[2]*e2[0] - e1[0]*e2[2];
                double cz = e1[0]*e2[1] - e1[1]*e2[0];
                faceAreaSum += 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
            }
        }
        vol(e) = mesh.GetElementVolume(e);
        dtscale(e) = vol(e) / (faceAreaSum / 2.0);
    }
    return dtscale;
}

double getJacobiGQ_RMin(const int order, const int basis_type) {
    auto mesh{ mfem::Mesh::MakeCartesian1D(1, 2.0) };
    mfem::DG_FECollection fec{ order, 1, basis_type };
    mfem::FiniteElementSpace fes{ &mesh, &fec };

    mfem::GridFunction nodes(&fes);
    mesh.GetNodes(nodes);

    return std::abs(nodes(0) - nodes(1));
}
std::vector<Source::Position> getVerticesCoordsForElem(const mfem::ParFiniteElementSpace& fes, const ElementId& e, const std::vector<Source::Position>& positions)
{
    mfem::Array<int> vertices;
    fes.GetElementVertices(e, vertices);
    std::vector<Source::Position> res(vertices.Size());
    for (auto v{ 0 }; v < vertices.Size(); v++) {
        res[v] = positions[vertices[v]];
    }
    return res;
}

double getSideLength(const Source::Position& va, const Source::Position& vb) {
    return std::sqrt((vb[0] - va[0]) * (vb[0] - va[0]) + (vb[1] - va[1]) * (vb[1] - va[1]));
}

double getElementPerimeter(const std::vector<Source::Position>& vertCoords)
{
    double res = 0.0;
    int n = vertCoords.size();
    for (int i = 0; i < n; i++) {
        res += getSideLength(vertCoords[i], vertCoords[(i + 1) % n]);
    }
    return res;
}

double calcMeshTimeStep(mfem::FiniteElementSpace& fes, double heuristic_divisor,
                        const mfem::Element::Type& special_elem_type, const int basis_type)
{
    Vector dtscale{ getTimeStepScale(*fes.GetMesh()) };
    double rmin{ getJacobiGQ_RMin(fes.FEColl()->GetOrder(), basis_type) };
    auto dt{ dtscale.Min() * rmin * 2.0 / 3.0 / physicalConstants::speedOfLight };
    dt *= 0.75; // Purely heuristic.
    if (checkIfElemTypeInMesh(*fes.GetMesh(), special_elem_type)) {
        return dt / heuristic_divisor;
    }
    return dt;
}

double estimateTimeStep(const Model& model, const SolverOptions& opts, const mfem::ParFiniteElementSpace& fes, const TimeDependentOperator* tdo)
{
    mfem::Mesh serialmesh = mfem::Mesh(model.getConstSerialMesh());
    mfem::DG_FECollection fec(fes.FEColl()->GetOrder(), serialmesh.Dimension(), opts.basis_type);
    mfem::FiniteElementSpace serialfes(&serialmesh, &fec);

    int dimension = model.getConstMesh().Dimension();

    // 1D dimension handling
    if (dimension == 1) {
        double maxTimeStep{ 0.0 };
        if (opts.evolution.order == 0) {
            maxTimeStep = getMinimumInterNodeDistance(serialfes) / physicalConstants::speedOfLight;
        }
        else {
            maxTimeStep = getMinimumInterNodeDistance(serialfes) / pow(double(fes.FEColl()->GetOrder()), 1.5) / physicalConstants::speedOfLight;
        }
        return opts.cfl * maxTimeStep;
    }

    // 2D and 3D common calculation
    double base_dt = 0.0;
    if (dimension == 2) {
        base_dt = calcMeshTimeStep(serialfes, 2.0, mfem::Element::Type::QUADRILATERAL, opts.basis_type) * opts.cfl;
    }
    else if (dimension == 3) {
        base_dt = calcMeshTimeStep(serialfes, 6.0, mfem::Element::Type::HEXAHEDRON, opts.basis_type) * opts.cfl;
    }
    else {
        throw std::runtime_error("Automatic Time Step Estimation not available for the set dimension.");
    }

    // Hesthaven-specific 3D scaling
    if (opts.evolution.op == EvolutionOperatorType::Hesthaven && dimension == 3) {
        auto maxFscaleVal{ 0.0 };
        const auto& evol = dynamic_cast<const HesthavenEvolution*>(tdo);
        for (auto e{ 0 }; e < serialmesh.GetNE(); e++) {
            const auto& fscaleMax{ evol->getHesthavenElement(e).fscale.maxCoeff() };
            if (fscaleMax > maxFscaleVal) {
                maxFscaleVal = fscaleMax;
            }
        }
        const auto& order = serialfes.FEColl()->GetOrder();
        return 1.0 * opts.cfl / (maxFscaleVal * order * order);
    }

    // Global-specific 3D scaling
    if (opts.evolution.op == EvolutionOperatorType::Global && dimension == 3) {
        return base_dt / 0.8; // 0.8 is purely heuristic, adjusted from Hesthaven ATS value.
    }

    return base_dt;
}

double Solver::calcAverageElementSizeInMesh()
{
    double res = 0.0;
    auto& mesh = this->model_.getMesh();

    for (int e = 0; e < mesh.GetNE(); e++)
    {
        res += mesh.GetElementSize(e);
    }

    return res / mesh.GetNE();
}

size_t getCurrentMemoryUsage() {
#ifdef SEMBA_DGTD_ENABLE_CUDA
    size_t free_bytes = 0;
    size_t total_bytes = 0;
    cudaError_t err = cudaMemGetInfo(&free_bytes, &total_bytes);
    if (err != cudaSuccess) {
        std::cerr << "CUDA memory query failed: " << cudaGetErrorString(err) << "\n";
        return 0;
    }
    return total_bytes - free_bytes; // bytes currently in use on GPU
#else
    std::ifstream statm("/proc/self/statm");
    long rss_pages = 0;
    statm >> rss_pages >> rss_pages;
    return static_cast<size_t>(rss_pages * sysconf(_SC_PAGESIZE));
#endif
}

size_t getPeakMemoryUsage() {
#ifdef SEMBA_DGTD_ENABLE_CUDA
    return getCurrentMemoryUsage(); // No peak tracking available; report current GPU usage.
#else
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.rfind("VmHWM:", 0) == 0) {
            long kb = 0;
            std::istringstream(line.substr(6)) >> kb;
            return static_cast<size_t>(kb * 1024);
        }
    }
    return 0;
#endif
}

void Solver::writeSimulationStatistics(const Time runtime){
    // Skip writing statistics for SGBC sub-solvers
    if (opts_.is_sgbc_solver) {
        return;
    }

    MPI_Comm comm = model_.getMesh().GetComm();
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::filesystem::path simExpPath(getSimulationCaseExportPath(model_.meshName_) + "/SimulationStats/");

    std::filesystem::create_directories(simExpPath);

    std::string path = (simExpPath / ("statistics_rank" + std::to_string(rank) + ".dat")).string();

    std::string existing;
    {
        std::ifstream previous(path);
        if (previous) {
            std::ostringstream buffered;
            buffered << previous.rdbuf();
            existing = buffered.str();
        }
    }
    const auto marker = existing.find("Simulation Run Time:");
    if (marker != std::string::npos) {
        const auto line_start = existing.rfind('\n', marker);
        existing.resize(line_start == std::string::npos ? 0 : line_start + 1);
    }

    std::ofstream myfile(path, std::ios::trunc);
    if (myfile.is_open()) {
        if (!existing.empty()) {
            myfile << existing;
            if (existing.back() != '\n') {
                myfile << '\n';
            }
        }
        std::size_t mem_baseline = 0;
        std::size_t mem_peak = 0;
        std::size_t mem_count = 0;
        double mem_sum = 0.0;
        snapshotTemporal(mem_baseline, mem_peak, mem_sum, mem_count);
        const std::size_t temporalMemAverageSampled = mem_count > 0
            ? static_cast<std::size_t>(mem_sum / static_cast<double>(mem_count))
            : 0;
        const std::size_t temporalIncrease = mem_peak > mem_baseline ? mem_peak - mem_baseline : 0;
        myfile << std::scientific << std::setprecision(5);
        myfile << "Simulation Run Time: " << runtime << " (s)\n";
        myfile << std::defaultfloat;
        myfile << "Final Time: " << (opts_.final_time / physicalConstants::speedOfLight_SI * 1e9) << " (ns)\n";
        myfile << "Time Step: " << (dt_ / physicalConstants::speedOfLight_SI * 1e9) << " (ns)\n";
        myfile << "Number of Mesh Elements: " << fes_.get()->GetNE() << "\n";
        myfile << "Average Element Size in Mesh: " << calcAverageElementSizeInMesh() << "\n";
        auto local_dofs = static_cast<std::int64_t>(fes_->GetNDofs());
        myfile << "Number of Local Degrees of Freedom: " << local_dofs << "\n";
        if (opts_.evolution.op == Global) {
            auto global = dynamic_cast<GlobalEvolution*>(evolTDO_.get());
            std::int64_t local_elems = static_cast<std::int64_t>(std::pow(global->getConstGlobalOperator().Size(), 2.0));
            myfile << "Operator Total Elements for Local Degrees of Freedom: " << local_elems << "\n";
            std::int64_t local_and_ghost_elems = static_cast<std::int64_t>(global->getConstGlobalOperator().Height()) * static_cast<std::int64_t>(global->getConstGlobalOperator().Width());
            myfile << "Operator Total Elements for Local and Ghost Degrees of Freedom: " << local_and_ghost_elems << "\n";
            std::int64_t non_zero_elems = static_cast<std::int64_t>(global->getConstGlobalOperator().NumNonZeroElems());
            myfile << "Number of Operator Non-Zero Elements: " << non_zero_elems << "\n";
        }
        myfile << "Temporal Evolution Baseline Memory (B): " << mem_baseline << "\n";
        myfile << "Temporal Evolution Peak Sampled Memory (B): " << mem_peak << "\n";
        myfile << "Temporal Evolution Average Sampled Memory (B): " << temporalMemAverageSampled << "\n";
        myfile << "Temporal Evolution Sample Count: " << mem_count << "\n";
        myfile << "Temporal Evolution Peak Increase (B): " << temporalIncrease << "\n";
        myfile << "Temporal Evolution Memory Consumption (B): " << temporalMemAverageSampled << "\n";
        myfile.close();
        const int fd = ::open(path.c_str(), O_RDONLY);
        if (fd >= 0) {
            ::fsync(fd);
            ::close(fd);
        }
    } else {
        std::cerr << "Rank " << rank << " failed to open file: " << path << "\n";
    }
}

#ifdef SHOW_TIMER_INFORMATION
void printSimulationInformation(const double time, const double dt, const double final_time)
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Information is updated every 30 seconds." << std::endl;
    std::cout << "Current Step: " + std::to_string(int(time / dt)) << std::endl;
    std::cout << "Steps Left  : " + std::to_string(int((final_time - time) / dt)) << std::endl;
    std::cout << std::endl;
    std::cout << "Final Time  : " + std::to_string(final_time / physicalConstants::speedOfLight_SI * 1e9) + " ns." << std::endl;
    std::cout << "Current Time: " + std::to_string(time / physicalConstants::speedOfLight_SI * 1e9) + " ns." << std::endl;
    std::cout << "Time Step   : " + std::to_string(dt / physicalConstants::speedOfLight_SI * 1e9) + " ns." << std::endl;
    std::cout << std::endl;
}
#endif

void loadSGBCNodalFieldsWithSolverFields(SGBCForcingFields& fields_in, const Fields<mfem::ParFiniteElementSpace,mfem::ParGridFunction>& global_fields, const SGBCGlobalNodeInfo& node_pairs)
{
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            fields_in.at(f).at(d).first  = global_fields.get(f,d)[node_pairs.g_el1];
            fields_in.at(f).at(d).second = global_fields.get(f,d)[node_pairs.g_el2];
        }
    }
}

void loadSolverFieldsWithSGBCValues(const SGBCForcingFields& fields_in, Fields<mfem::ParFiniteElementSpace,mfem::ParGridFunction>& global_fields, const SGBCGlobalNodeInfo& node_pairs)
{
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            global_fields.get(f,d)[node_pairs.g_el1] += fields_in.at(f).at(d).first;
            global_fields.get(f,d)[node_pairs.g_el2] += fields_in.at(f).at(d).second;
        }
    }
}

namespace {

std::atomic<int> g_checkpoint_stop{0};

void onCheckpointSignal(int)
{
    g_checkpoint_stop.store(1, std::memory_order_relaxed);
}

class CheckpointSignalGuard {
public:
    CheckpointSignalGuard()
    {
        g_checkpoint_stop.store(0, std::memory_order_relaxed);
        struct sigaction action {};
        action.sa_handler = onCheckpointSignal;
        sigemptyset(&action.sa_mask);
        action.sa_flags = SA_RESTART;
        sigaction(SIGINT, &action, &old_int_);
        sigaction(SIGTERM, &action, &old_term_);
    }

    ~CheckpointSignalGuard()
    {
        sigaction(SIGINT, &old_int_, nullptr);
        sigaction(SIGTERM, &old_term_, nullptr);
    }

    CheckpointSignalGuard(const CheckpointSignalGuard&) = delete;
    CheckpointSignalGuard& operator=(const CheckpointSignalGuard&) = delete;

private:
    struct sigaction old_int_ {};
    struct sigaction old_term_ {};
};

struct SgbcKey {
    int tag = 0;
    int node_a = 0;
    int node_b = 0;

    bool operator<(const SgbcKey& other) const
    {
        if (tag != other.tag) return tag < other.tag;
        if (node_a != other.node_a) return node_a < other.node_a;
        return node_b < other.node_b;
    }
};

}  // namespace

void Solver::enableCheckpoint(const std::string& json_path, const std::string& mesh_path)
{
    checkpoint_json_path_ = json_path;
    checkpoint_mesh_path_ = mesh_path;
}

double Solver::accumulatedRunSeconds() const
{
    double segment = 0.0;
    if (run_clock_started_) {
        segment = std::chrono::duration<double>(std::chrono::steady_clock::now() - run_start_).count();
    }
    return checkpoint_elapsed_run_ + segment;
}

void Solver::snapshotTemporal(std::size_t& baseline, std::size_t& peak, double& sum, std::size_t& count) const
{
    baseline = checkpoint_temporal_count_ > 0 ? checkpoint_temporal_baseline_ : temporalMemBaseline_;
    peak = std::max(checkpoint_temporal_peak_, temporalMemPeakSampled_);
    sum = checkpoint_temporal_sum_ + static_cast<double>(temporalMemSumSampled_);
    count = checkpoint_temporal_count_ + temporalMemSampleCount_;
}

void Solver::foldTemporalAndClock(double elapsed_seconds, std::chrono::steady_clock::time_point stamp)
{
    std::size_t baseline = 0;
    std::size_t peak = 0;
    std::size_t count = 0;
    double sum = 0.0;
    snapshotTemporal(baseline, peak, sum, count);
    checkpoint_temporal_baseline_ = baseline;
    checkpoint_temporal_peak_ = peak;
    checkpoint_temporal_sum_ = sum;
    checkpoint_temporal_count_ = count;
    temporalMemBaseline_ = 0;
    temporalMemPeakSampled_ = 0;
    temporalMemSumSampled_ = 0.0L;
    temporalMemSampleCount_ = 0;
    checkpoint_elapsed_run_ = elapsed_seconds;
    run_start_ = stamp;
    run_clock_started_ = true;
}

void Solver::purgeCheckpoints()
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        try {
            const auto root = checkpointRoot(model_.meshName_);
            if (std::filesystem::exists(root)) {
                std::filesystem::remove_all(root);
            }
        } catch (const std::exception& ex) {
            std::cerr << "Finished run, but checkpoints could not be removed: " << ex.what() << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

bool Solver::writeCheckpoint()
{
    if (checkpoint_json_path_.empty() || checkpoint_mesh_path_.empty()) {
        return false;
    }

    flushProbeFiles();

    const auto stamp = std::chrono::steady_clock::now();
    const double local_elapsed = run_clock_started_
        ? checkpoint_elapsed_run_ + std::chrono::duration<double>(stamp - run_start_).count()
        : checkpoint_elapsed_run_;
    double elapsed = local_elapsed;
    if (Mpi::WorldSize() > 1) {
        MPI_Allreduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    std::size_t mem_baseline = 0;
    std::size_t mem_peak = 0;
    std::size_t mem_count = 0;
    double mem_sum = 0.0;
    snapshotTemporal(mem_baseline, mem_peak, mem_sum, mem_count);

    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    std::vector<unsigned long long> baselines(static_cast<std::size_t>(nprocs));
    std::vector<unsigned long long> peaks(static_cast<std::size_t>(nprocs));
    std::vector<unsigned long long> counts(static_cast<std::size_t>(nprocs));
    std::vector<double> sums(static_cast<std::size_t>(nprocs));
    const unsigned long long local_baseline = mem_baseline;
    const unsigned long long local_peak = mem_peak;
    const unsigned long long local_count = mem_count;
    MPI_Gather(&local_baseline, 1, MPI_UNSIGNED_LONG_LONG, baselines.data(), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_peak, 1, MPI_UNSIGNED_LONG_LONG, peaks.data(), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_count, 1, MPI_UNSIGNED_LONG_LONG, counts.data(), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Gather(&mem_sum, 1, MPI_DOUBLE, sums.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    CheckpointManifest meta;
    meta.time = time_;
    meta.dt = dt_;
    meta.final_time = opts_.final_time;
    meta.order = opts_.evolution.order;
    meta.cycle = probesManager_.cycle();
    meta.next_checkpoint_mark = checkpoint_next_mark_;
    meta.elapsed_run_seconds = elapsed;
    meta.json_sha256 = checkpoint_json_sha256_;
    meta.mesh_sha256 = checkpoint_mesh_sha256_;
    meta.exporters = probesManager_.captureExporterCursors();
    meta.mor = probesManager_.captureMorCursors();
    if (rank == 0) {
        if (checkpoint_json_sha256_.empty()) {
            checkpoint_json_sha256_ = sha256File(checkpoint_json_path_);
            checkpoint_mesh_sha256_ = sha256File(checkpoint_mesh_path_);
            meta.json_sha256 = checkpoint_json_sha256_;
            meta.mesh_sha256 = checkpoint_mesh_sha256_;
        }
        meta.temporal_mem_baseline.assign(baselines.begin(), baselines.end());
        meta.temporal_mem_peak.assign(peaks.begin(), peaks.end());
        meta.temporal_mem_count.assign(counts.begin(), counts.end());
        meta.temporal_mem_sum = std::move(sums);
    }

    auto& state = fields_.allDOFs();
    const double* state_ptr = state.HostRead();

    std::vector<SgbcRecord> sgbc;
    if (globalEvol_cache_ != nullptr) {
        for (const auto& [tag, states] : globalEvol_cache_->sgbcStates()) {
            for (const auto& face : states) {
                SgbcRecord record;
                record.tag = tag;
                record.node_a = face.global_pair.first;
                record.node_b = face.global_pair.second;
                const int n = face.fields_state.Size();
                if (n > 0) {
                    const double* values = face.fields_state.HostRead();
                    record.values.assign(values, values + n);
                }
                sgbc.push_back(std::move(record));
            }
        }
    }

    const auto& partition = model_.elementPartition();
    try {
        writeCheckpointCollective(
            model_.meshName_,
            checkpoint_json_path_,
            checkpoint_mesh_path_,
            meta,
            partition.GetData(),
            partition.Size(),
            state_ptr,
            state.Size(),
            sgbc);
    } catch (const std::exception& ex) {
        checkpoint_hold_ = true;
        checkpoint_fail_at_ = std::chrono::steady_clock::now();
        if (rank == 0) {
            std::cerr << "Checkpoint failed (" << ex.what()
                      << "). The run will continue; the last complete checkpoint is unchanged."
                      << std::endl;
        }
        return false;
    }
    checkpoint_hold_ = false;
    foldTemporalAndClock(elapsed, stamp);
    return true;
}

void Solver::restoreCheckpoint(const std::string& directory)
{
    const MPI_Comm comm = model_.getMesh().GetComm();
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::string error;
    int failed = 0;
    try {
        const auto manifest = readManifest(std::filesystem::path(directory) / "manifest.json");
        if (manifest.format_version != kCheckpointFormatVersion) {
            throw std::runtime_error(
                "Checkpoint format version " + std::to_string(manifest.format_version)
                + " is not supported (this build reads version "
                + std::to_string(kCheckpointFormatVersion) + ").");
        }
        if (manifest.world_size != Mpi::WorldSize()) {
            throw std::runtime_error(
                "Checkpoint was written with " + std::to_string(manifest.world_size)
                + " ranks; this job has " + std::to_string(Mpi::WorldSize()) + ".");
        }
        const auto state_path = std::filesystem::path(directory) / ("state.rank" + std::to_string(rank) + ".bin");
        if (!manifest.state_sha256.empty()) {
            if (rank < 0 || rank >= static_cast<int>(manifest.state_sha256.size())
                || sha256File(state_path) != manifest.state_sha256[static_cast<std::size_t>(rank)]) {
                throw std::runtime_error(
                    "Checkpoint state file for rank " + std::to_string(rank) + " failed its checksum.");
            }
        }
        if (rank < 0 || rank >= static_cast<int>(manifest.state_sizes.size())
            || manifest.state_sizes[static_cast<std::size_t>(rank)] != fields_.allDOFs().Size()) {
            const int saved = (rank >= 0 && rank < static_cast<int>(manifest.state_sizes.size()))
                ? manifest.state_sizes[static_cast<std::size_t>(rank)] : -1;
            throw std::runtime_error(
                "Checkpoint state length on rank " + std::to_string(rank)
                + " is " + std::to_string(saved)
                + "; this mesh partition has " + std::to_string(fields_.allDOFs().Size()) + ".");
        }

        std::vector<SgbcRecord> saved_sgbc;
        readRankState(
            state_path,
            fields_.allDOFs().HostWrite(),
            fields_.allDOFs().Size(),
            saved_sgbc);

        int live_sgbc = 0;
        if (globalEvol_cache_ != nullptr) {
            for (const auto& [tag, states] : globalEvol_cache_->sgbcStates()) {
                live_sgbc += static_cast<int>(states.size());
                (void)tag;
            }
        }
        const int saved_count = (rank < static_cast<int>(manifest.sgbc_counts.size()))
            ? manifest.sgbc_counts[static_cast<std::size_t>(rank)] : -1;
        if (saved_count != live_sgbc || static_cast<int>(saved_sgbc.size()) != live_sgbc) {
            throw std::runtime_error(
                "Checkpoint SGBC state count on rank " + std::to_string(rank)
                + " is " + std::to_string(saved_count)
                + "; this run has " + std::to_string(live_sgbc) + ".");
        }
        if (live_sgbc > 0) {
            std::map<SgbcKey, std::size_t> index;
            for (std::size_t i = 0; i < saved_sgbc.size(); ++i) {
                SgbcKey key;
                key.tag = saved_sgbc[i].tag;
                key.node_a = saved_sgbc[i].node_a;
                key.node_b = saved_sgbc[i].node_b;
                index.emplace(key, i);
            }
            int matched = 0;
            for (auto& [tag, states] : globalEvol_cache_->sgbcStates()) {
                for (auto& face : states) {
                    SgbcKey key;
                    key.tag = tag;
                    key.node_a = face.global_pair.first;
                    key.node_b = face.global_pair.second;
                    const auto it = index.find(key);
                    if (it == index.end()) {
                        throw std::runtime_error("Checkpoint is missing an SGBC face state.");
                    }
                    const auto& values = saved_sgbc[it->second].values;
                    if (static_cast<int>(values.size()) != face.fields_state.Size()) {
                        throw std::runtime_error("Checkpoint SGBC face state has the wrong length.");
                    }
                    if (!values.empty()) {
                        std::memcpy(
                            face.fields_state.HostWrite(),
                            values.data(),
                            values.size() * sizeof(double));
                    }
                    ++matched;
                }
            }
            if (matched != live_sgbc) {
                throw std::runtime_error("Checkpoint SGBC states do not match this mesh.");
            }
        }

        time_ = manifest.time;
        dt_ = manifest.dt;
        if (evolTDO_) {
            evolTDO_->SetTime(time_);
        }
        checkpoint_next_mark_ = std::max(1, manifest.next_checkpoint_mark);
        checkpoint_elapsed_run_ = manifest.elapsed_run_seconds;
        checkpoint_json_sha256_ = manifest.json_sha256;
        checkpoint_mesh_sha256_ = manifest.mesh_sha256;
        if (rank >= 0
            && rank < static_cast<int>(manifest.temporal_mem_count.size())
            && rank < static_cast<int>(manifest.temporal_mem_baseline.size())
            && rank < static_cast<int>(manifest.temporal_mem_peak.size())
            && rank < static_cast<int>(manifest.temporal_mem_sum.size())) {
            const auto idx = static_cast<std::size_t>(rank);
            checkpoint_temporal_baseline_ = static_cast<std::size_t>(manifest.temporal_mem_baseline[idx]);
            checkpoint_temporal_peak_ = static_cast<std::size_t>(manifest.temporal_mem_peak[idx]);
            checkpoint_temporal_sum_ = manifest.temporal_mem_sum[idx];
            checkpoint_temporal_count_ = static_cast<std::size_t>(manifest.temporal_mem_count[idx]);
        }
        probesManager_.recalculateExportSteps(dt_);
        probesManager_.restoreCheckpointCursors(manifest.cycle, manifest.exporters, manifest.mor);
        probesManager_.trimProbeOutput(time_);
        if (rank == 0) {
            std::cout << "Resuming from " << directory << " at t = " << time_ << std::endl;
        }
    } catch (const std::exception& ex) {
        failed = 1;
        error = ex.what();
    }

    int global_failed = 0;
    MPI_Allreduce(&failed, &global_failed, 1, MPI_INT, MPI_MAX, comm);
    if (global_failed) {
        if (!error.empty()) {
            throw std::runtime_error(error);
        }
        throw std::runtime_error("Checkpoint load failed on another rank.");
    }
}

void Solver::run()
{

    run_start_ = std::chrono::steady_clock::now();
    run_clock_started_ = true;
    temporalMemBaseline_ = getCurrentMemoryUsage();
    temporalMemPeakSampled_ = temporalMemBaseline_;
    temporalMemSumSampled_ = static_cast<long double>(temporalMemBaseline_);
    temporalMemSampleCount_ = 1;

#ifdef SHOW_TIMER_INFORMATION
    auto lastPrintTime{ std::chrono::steady_clock::now() };
    if (Mpi::WorldRank() == 0){
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "-------------SOLVER RUN INFORMATION-------------" << std::endl;
        printSimulationInformation(time_, dt_, opts_.final_time);
    }
#endif

    const bool checkpointing = !checkpoint_json_path_.empty();
    bool stopped_by_signal = false;
    std::optional<CheckpointSignalGuard> checkpoint_signals;
    if (checkpointing) {
        checkpoint_signals.emplace();
    }

    while (time_ <= opts_.final_time - 1e-8*dt_) {
        step();
        sampleTemporalMemory();

        // Stability check — every step (getNorml2 is O(N), negligible cost).
        // NaN > threshold is always false in IEEE 754, so we must use isfinite.
        {
            double localNorm = this->fields_.getNorml2();
            int localUnstable = (!std::isfinite(localNorm) || localNorm > 1e20) ? 1 : 0;
            int globalUnstable = 0;
            MPI_Allreduce(&localUnstable, &globalUnstable, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (globalUnstable && Mpi::WorldRank() == 0) {
                std::cout << "========================================================================" << std::endl;
                std::cout << "WARNING: Simulation is potentially unstable (field norm = "
                          << localNorm << ")." << std::endl;
                std::cout << "  Time: " << time_ << " / " << opts_.final_time
                          << ",  dt: " << dt_ << std::endl;
                std::cout << "  Verify your setup and consider lowering the time step." << std::endl;
                std::cout << "========================================================================" << std::endl;
            }
        }

#ifdef SHOW_TIMER_INFORMATION
        auto currentTime = std::chrono::steady_clock::now();
        int local_due = (std::chrono::duration_cast<std::chrono::seconds>
            (currentTime - lastPrintTime).count() >= 30.0) ? 1 : 0;
        int global_due = local_due;
        if (Mpi::WorldSize() > 1) {
            MPI_Allreduce(&local_due, &global_due, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (global_due)
        {
            if (Mpi::WorldRank() == 0){
                printSimulationInformation(time_, dt_, opts_.final_time);
            }
            const int P = Mpi::WorldSize();
            int local_gather = (!opts_.is_sgbc_solver && stepTimingStats_.step_count > 0) ? 1 : 0;
            int global_gather = local_gather;
            if (P > 1) {
                MPI_Allreduce(&local_gather, &global_gather, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (global_gather) {
                double local_avg[5] = {0, 0, 0, 0, 0};
                if (local_gather) {
                    const double n = static_cast<double>(stepTimingStats_.step_count);
                    local_avg[0] = stepTimingStats_.step_ms / n;
                    local_avg[1] = stepTimingStats_.ode_ms / n;
                    local_avg[2] = stepTimingStats_.sgbc_finalize_ms / n;
                    local_avg[3] = stepTimingStats_.probe_sync_ms / n;
                    local_avg[4] = stepTimingStats_.probe_update_ms / n;
                }
                std::vector<double> all_avg(5 * P);
                if (P > 1) {
                    MPI_Gather(local_avg, 5, MPI_DOUBLE,
                               all_avg.data(), 5, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
                } else {
                    std::copy(local_avg, local_avg + 5, all_avg.data());
                }

                if (Mpi::WorldRank() == 0) {
                    std::cout << "[Step timing] avg of " << stepTimingStats_.step_count
                              << " steps, ms/step\n";
                    std::cout << "  Rank  | total     ode  sgbc_fin  probe_sync  probe_update\n";
                    std::cout << "  ------+---------------------------------------------------\n";
                    for (int r = 0; r < P; ++r) {
                        std::cout << std::setw(6) << r << " | "
                                  << std::setw(6) << all_avg[r*5 + 0] << " "
                                  << std::setw(7) << all_avg[r*5 + 1] << " "
                                  << std::setw(9) << all_avg[r*5 + 2] << " "
                                  << std::setw(11) << all_avg[r*5 + 3] << " "
                                  << std::setw(12) << all_avg[r*5 + 4] << "\n";
                        }
                }
                stepTimingStats_ = StepTimingStats{};
                probesManager_.printTimingSummaryAndReset();
            }
            lastPrintTime = currentTime;
        }
#endif
        if (checkpointing) {
            int local_stop = g_checkpoint_stop.load(std::memory_order_relaxed);
            int global_stop = local_stop;
            if (Mpi::WorldSize() > 1) {
                MPI_Allreduce(&local_stop, &global_stop, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }

            int next_mark = checkpoint_next_mark_;
            int due = 0;
            if (Mpi::WorldRank() == 0) {
                const CheckpointSchedule schedule = checkpointSchedule(
                    time_, dt_, opts_.final_time, opts_.checkpoint_percent, next_mark);
                due = schedule.due ? 1 : 0;
                next_mark = schedule.next_mark;
                if (due && checkpoint_hold_) {
                    const double since = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - checkpoint_fail_at_).count();
                    if (since < 60.0) {
                        due = 0;
                    }
                }
            }
            if (Mpi::WorldSize() > 1) {
                MPI_Bcast(&due, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&next_mark, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
            if (global_stop && Mpi::WorldRank() == 0) {
                std::cout << "Received stop signal. Writing checkpoint and exiting." << std::endl;
            }
            if (global_stop || due) {
                const int saved_mark = checkpoint_next_mark_;
                if (due) {
                    checkpoint_next_mark_ = next_mark;
                }
                if (!writeCheckpoint()) {
                    checkpoint_next_mark_ = saved_mark;
                }
            }
            if (global_stop) {
                stopped_by_signal = true;
                break;
            }
        }
    }

    writeSimulationStatistics(accumulatedRunSeconds());
    if (checkpointing && !stopped_by_signal) {
        purgeCheckpoints();
    }

}

void Solver::step(bool update_probes)
{
#ifdef SHOW_TIMER_INFORMATION
    using clock = std::chrono::steady_clock;
    auto step_start = clock::now();
#endif
    double truedt{ std::min(dt_, opts_.final_time - time_) };

    // Monolithic IMEX: checkpoint SGBC before RK4, advance inside Mult(), finalize after
    if (globalEvol_cache_ && globalEvol_cache_->hasSGBC()) {
        globalEvol_cache_->commitSGBCCheckpoint(time_, truedt, fields_);
    }

#ifdef SHOW_TIMER_INFORMATION
    auto ode_start = clock::now();
#endif
    odeSolver_->Step(fields_.allDOFs(), time_, truedt);
#ifdef SEMBA_DGTD_ENABLE_CUDA
    // Make ODE timing honest: RK4 Mult kernels are async on CUDA.
    // Do not use MFEM_STREAM_SYNC here — it is empty unless compiled with nvcc.
    if (mfem::Device::Allows(mfem::Backend::CUDA)) {
        cudaStreamSynchronize(0);
    }
#endif
#ifdef SHOW_TIMER_INFORMATION
    auto after_ode = clock::now();
#endif

    if (globalEvol_cache_ && globalEvol_cache_->hasSGBC()) {
        globalEvol_cache_->finalizeSGBCStep(fields_);
    }
#ifdef SHOW_TIMER_INFORMATION
    auto after_finalize = clock::now();
#endif

    if (update_probes) {
#ifdef SHOW_TIMER_INFORMATION
        auto probe_sync_start = clock::now();
#endif
        // Full-volume HostRead / ExchangeFaceNbrData removed from the hot path.
        // Consumers pull only what they need:
        //   - RCS: device gather of surface DOFs + small HostRead
        //   - Exporter / point / field: HostRead inside their update paths
#ifdef SHOW_TIMER_INFORMATION
        auto after_probe_sync = clock::now();
#endif
        probesManager_.updateProbes(time_);
#ifdef SHOW_TIMER_INFORMATION
        auto after_probe_update = clock::now();
#endif
#ifdef SHOW_TIMER_INFORMATION
        stepTimingStats_.step_ms += std::chrono::duration<double, std::milli>(after_probe_update - step_start).count();
        stepTimingStats_.ode_ms += std::chrono::duration<double, std::milli>(after_ode - ode_start).count();
        stepTimingStats_.sgbc_finalize_ms += std::chrono::duration<double, std::milli>(after_finalize - after_ode).count();
        stepTimingStats_.probe_sync_ms += std::chrono::duration<double, std::milli>(after_probe_sync - probe_sync_start).count();
        stepTimingStats_.probe_update_ms += std::chrono::duration<double, std::milli>(after_probe_update - after_probe_sync).count();
        stepTimingStats_.step_count++;
#endif
    } else {
#ifdef SHOW_TIMER_INFORMATION
        auto step_end = clock::now();
        stepTimingStats_.step_ms += std::chrono::duration<double, std::milli>(step_end - step_start).count();
        stepTimingStats_.ode_ms += std::chrono::duration<double, std::milli>(after_ode - ode_start).count();
        stepTimingStats_.sgbc_finalize_ms += std::chrono::duration<double, std::milli>(after_finalize - after_ode).count();
        stepTimingStats_.step_count++;
#endif
    }
}


GeomTagToBoundary Solver::assignAttToBdrByDimForSpectral(mfem::ParMesh& submesh)
{
    switch (submesh.Dimension()) {
    case 1:
        return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA} };
    case 2:
        switch (submesh.GetElementType(0)) {
        case mfem::Element::TRIANGLE:
            return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA } };
        case mfem::Element::QUADRILATERAL:
            return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA} };
        default:
            throw std::runtime_error("Incorrect element type for 2D spectral AttToBdr assignation.");
        }
    case 3:
        switch (submesh.GetElementType(0)) {
        case mfem::Element::TETRAHEDRON:
            return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA} };
        case mfem::Element::HEXAHEDRON:
            return GeomTagToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA}, {5, BdrCond::SMA }, {6, BdrCond::SMA} };
        default:
            throw std::runtime_error("Incorrect element type for 3D spectral AttToBdr assignation.");
        }
    default:
        throw std::runtime_error("Dimension is incorrect for spectral AttToBdr assignation.");
    }

}

Eigen::SparseMatrix<double> Solver::assembleSubmeshedSpectralOperatorMatrix(mfem::ParMesh& submesh, const mfem::FiniteElementCollection& fec, const EvolutionOptions& opts)
{
    Model submodel(submesh, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo(assignAttToBdrByDimForSpectral(submesh), GeomTagToInteriorBoundary{}));
    mfem::ParFiniteElementSpace subfes(&submesh, &fec);
    Eigen::SparseMatrix<double> local;
    auto numberOfFieldComponents = 2;
    auto numberofMaxDimensions = 3;
    local.resize(numberOfFieldComponents * numberofMaxDimensions * subfes.GetNDofs(), 
        numberOfFieldComponents * numberofMaxDimensions * subfes.GetNDofs());

    EvolutionOptions localopts(opts);
    ProblemDescription pd(submodel, probesManager_.probes, sourcesManager_.sources, localopts);
    DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, subfes);
    for (int x = X; x <= Z; x++) {
        int y = (x + 1) % 3;
        int z = (x + 2) % 3;

        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildDerivativeSubOperator<mfem::ParBilinearForm>(y)->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,z }, -1.0); // MS
        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildDerivativeSubOperator<mfem::ParBilinearForm>(z)->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,y });
        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildDerivativeSubOperator<mfem::ParBilinearForm>(y)->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,z });
        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildDerivativeSubOperator<mfem::ParBilinearForm>(z)->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,y }, -1.0);

        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildOneNormalSubOperator<mfem::ParBilinearForm>(E, { y })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,z }); // MFN
        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildOneNormalSubOperator<mfem::ParBilinearForm>(E, { z })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,y }, -1.0);
        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildOneNormalSubOperator<mfem::ParBilinearForm>(H, { y })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,z }, -1.0);
        allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildOneNormalSubOperator<mfem::ParBilinearForm>(H, { z })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,y });

        if (opts.alpha > 0.0) {

            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildZeroNormalSubOperator<mfem::ParBilinearForm>(H)->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { x }, -1.0); // MP
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildZeroNormalSubOperator<mfem::ParBilinearForm>(E)->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { x }, -1.0);
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildTwoNormalSubOperator<mfem::ParBilinearForm>(H, { X, x })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { X,x }); //MPNN
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildTwoNormalSubOperator<mfem::ParBilinearForm>(H, { Y, x })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { Y,x });
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(H)->SpMat(), dgops.buildTwoNormalSubOperator<mfem::ParBilinearForm>(H, { Z, x })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { Z,x });
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildTwoNormalSubOperator<mfem::ParBilinearForm>(E, { X, x })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { X,x });
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildTwoNormalSubOperator<mfem::ParBilinearForm>(E, { Y, x })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { Y,x });
            allocateDenseInEigen(buildByMult<mfem::ParFiniteElementSpace,mfem::ParBilinearForm>(dgops.buildInverseMassMatrixSubOperator<mfem::ParBilinearForm>(E)->SpMat(), dgops.buildTwoNormalSubOperator<mfem::ParBilinearForm>(E, { Z, x })->SpMat(), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { Z,x });

        }

    }
    return local;
}

double Solver::findMaxEigenvalueModulus(const Eigen::VectorXcd& eigvals)
{
    auto res{ 0.0 };
    for (int i = 0; i < eigvals.size(); ++i) {
        auto modulus{ sqrt(pow(eigvals[i].real(),2.0) + pow(eigvals[i].imag(),2.0)) };
        if (modulus <= 1.0 && modulus >= res) {
            res = modulus;
        }
    }
    return res;
}

void reassembleSpectralBdrForSubmesh(mfem::ParSubMesh* submesh) 
{
    switch (submesh->GetElementType(0)) {
    case mfem::Element::SEGMENT:
        for (int i = 0; i < submesh->GetParentVertexIDMap().Size(); ++i) {
            submesh->AddBdrPoint(i, i + 1);
        }
        submesh->FinalizeMesh();
        break;
    case mfem::Element::TRIANGLE:
        for (int i = 0; i < submesh->GetNBE(); ++i) {
            submesh->SetBdrAttribute(i, i + 1);
        }
        submesh->FinalizeMesh();
        break;
    case mfem::Element::QUADRILATERAL:
        for (int i = 0; i < submesh->GetNBE(); ++i) {
            submesh->SetBdrAttribute(i, i + 1);
        }
        submesh->FinalizeMesh();
        break;
    case mfem::Element::TETRAHEDRON:
        for (int i = 0; i < submesh->GetNBE(); ++i) {
            submesh->SetBdrAttribute(i, i + 1);
        }
        submesh->FinalizeMesh();
        break;
    case mfem::Element::HEXAHEDRON:
        for (int i = 0; i < submesh->GetNBE(); ++i) {
            submesh->SetBdrAttribute(i, i + 1);
        }
        submesh->FinalizeMesh();
        break;
    default:
        throw std::runtime_error("Incorrect element type for Bdr Spectral assignation.");
    }
}

void Solver::evaluateStabilityByEigenvalueEvolutionFunction(
    Eigen::VectorXcd& eigenvals, 
    MaxwellEvolution& maxwellEvol)
{
    auto real { toMFEMVector(eigenvals.real()) };
    auto realPre = real;
    auto imag { toMFEMVector(eigenvals.imag()) };
    auto imagPre = imag;
    auto time { 0.0 };
    maxwellEvol.SetTime(time);
    odeSolver_->Init(maxwellEvol);
    odeSolver_->Step(real, time, opts_.time_step);
    time = 0.0;
    maxwellEvol.SetTime(time);
    odeSolver_->Init(maxwellEvol);
    odeSolver_->Step(imag, time, opts_.time_step);
    
    for (int i = 0; i < real.Size(); ++i) {
        
        auto modPre{ sqrt(pow(realPre[i],2.0) + pow(imagPre[i],2.0)) };
        auto mod   { sqrt(pow(real[i]   ,2.0) + pow(imag[i]   ,2.0)) };

        if (modPre != 0.0) {
            if (mod / modPre > 1.0) {
                throw std::runtime_error("The coefficient between the modulus of a time evolved eigenvalue and its original value is higher than 1.0 - RK4 instability.");
            }
        }
    }
}

void Solver::performSpectralAnalysis(const mfem::ParFiniteElementSpace& fes, Model& model, const EvolutionOptions& opts)
{
    if (Mpi::WorldSize() > 1) {
        throw std::runtime_error(
            "Spectral analysis calls ParSubMesh once per local element and is serial only. Run with one MPI rank.");
    }
    mfem::Array<int> domainAtts(1);
    domainAtts[0] = 501;
    auto mesh{ model.getConstMesh() };
    auto meshCopy{ mesh };

    for (int elem = 0; elem < meshCopy.GetNE(); ++elem) {

        auto preAtt(meshCopy.GetAttribute(elem));
        meshCopy.SetAttribute(elem, domainAtts[0]);
        auto submesh{ mfem::ParSubMesh::CreateFromDomain(meshCopy,domainAtts) };
        meshCopy.SetAttribute(elem, preAtt);
        submesh.SetAttribute(0, preAtt);

        reassembleSpectralBdrForSubmesh(&submesh);

        auto eigenvals{ 
            assembleSubmeshedSpectralOperatorMatrix(submesh, *fes.FEColl(), opts).toDense().eigenvalues() 
        };
        mfem::ParFiniteElementSpace submeshFES{ &submesh, fes.FEColl() };
        Model model{ submesh,
            GeomTagToMaterialInfo{},
            GeomTagToBoundaryInfo(assignAttToBdrByDimForSpectral(submesh),GeomTagToInteriorBoundary{})
        };
        SourcesManager srcs{ Sources(), submeshFES, fields_ };
        ProblemDescription pd(model, probesManager_.probes, sourcesManager_.sources, opts_.evolution);
        MaxwellEvolution evol(pd, submeshFES, sourcesManager_);
        evaluateStabilityByEigenvalueEvolutionFunction(eigenvals, evol);
    }
}


}