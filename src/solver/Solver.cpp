#include "Solver.h"
#include "SolverExtension.h"
#include <filesystem>
#include <fstream>
#include <unistd.h>

namespace maxwell {

Solver::~Solver() = default; 

std::unique_ptr<mfem::ParFiniteElementSpace> buildFiniteElementSpace(mfem::ParMesh* m, mfem::FiniteElementCollection* fec)
{
    MPI_Barrier(MPI_COMM_WORLD);
    auto fes = std::make_unique<mfem::ParFiniteElementSpace>(m, fec);
    fes->ExchangeFaceNbrData();
    MPI_Barrier(MPI_COMM_WORLD);
    return fes;
}

std::unique_ptr<mfem::TimeDependentOperator> Solver::assignEvolutionOperator()
{
    if (opts_.evolution.op == EvolutionOperatorType::Hesthaven) {
        return std::make_unique<HesthavenEvolution>(*fes_, model_, sourcesManager_, opts_.evolution);
    }
    else if (opts_.evolution.op == EvolutionOperatorType::Global) {
        auto global_evol = std::make_unique<GlobalEvolution>(*fes_, model_, sourcesManager_, opts_.evolution);
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

        case ode_type::SDIRK23:  // L-stable flavor (good with PML/loss)
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
    fields_{ *fes_ },
    sourcesManager_{ sources, *fes_, fields_ },
    probesManager_ { probes , *fes_, fields_, opts_ },
    time_{0.0}
{
    auto initStartTime = std::chrono::steady_clock::now();

    std::filesystem::path simExpPath("Exports/" + getRunModeTag() + "/" + model_.meshName_ + "/SimulationStats/");

    if (Mpi::WorldRank() == 0){
        if (std::filesystem::exists(simExpPath)) {
            std::filesystem::remove_all(simExpPath);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::filesystem::create_directories(simExpPath);

    MPI_Barrier(MPI_COMM_WORLD);


    if (Mpi::WorldRank() == 0){
        checkOptionsAreValid(opts_);
    }

    if (opts_.evolution.spectral == true) {
        performSpectralAnalysis(*fes_.get(), model_, opts_.evolution);
    }

    evolTDO_ = assignEvolutionOperator();
    evolTDO_->SetTime(time_);

    if (opts_.time_step == 0.0) {
        dt_ = estimateTimeStep(model_, opts_, *fes_, evolTDO_.get());
    }
    else {
        dt_ = opts_.time_step;
    }

    assignODESolver();
    odeSolver_->Init(*evolTDO_);

    probesManager_.setCaseName(model_.meshName_);
    probesManager_.initPointFieldProbeExport();
    probesManager_.updateProbes(time_);

    auto initEndTime = std::chrono::steady_clock::now();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string path = (simExpPath / ("statistics_rank" + std::to_string(rank) + ".dat")).string();

    std::ofstream myfile(path, std::ios::app);
    if (myfile.is_open()) {
        auto runtime = std::chrono::duration<double>(initEndTime - initStartTime).count();
        myfile << std::scientific << std::setprecision(5);
        myfile << "Initialization Time: " << runtime << " (s)\n";
        myfile << std::defaultfloat;
        myfile.close();
    } else {
        std::cerr << "Rank " << rank << " failed to open file: " << path << "\n";
    }
}

void Solver::checkOptionsAreValid(const SolverOptions& opts) const
{
    
    if ((opts.evolution.order < 0) ||
        (opts.final_time < 0)) {
        throw std::runtime_error("Incorrect parameters in Options");
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
    mfem::Vector vol(mesh.GetNE()), dtscale(mesh.GetNE());
    for (int e = 0; e < mesh.GetNE(); ++e) {
        auto el{ mesh.GetElement(e) };
        mfem::Vector areasum(mesh.GetNumFaces());
        areasum = 0.0;
        for (int f = 0; f < mesh.GetElement(e)->GetNEdges(); ++f) {
            mfem::ElementTransformation* T{ mesh.GetFaceTransformation(f)};
            const mfem::IntegrationRule& ir = IntRules.Get(T->GetGeometryType(), T->OrderJ());
            for (int p = 0; p < ir.GetNPoints(); p++)
            {
                const mfem::IntegrationPoint& ip = ir.IntPoint(p);
                areasum(e) += ip.weight * T->Weight();
            }
        }
        vol(e) = mesh.GetElementVolume(e);
        dtscale(e) = vol(e) / (areasum(e) / 2.0);
    }
    return dtscale;
}

double getJacobiGQ_RMin(const int order) {
    auto mesh{ mfem::Mesh::MakeCartesian1D(1, 2.0) };
    mfem::DG_FECollection fec{ order,1, mfem::BasisType::GaussLobatto };
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
                        const mfem::Element::Type& special_elem_type)
{
    Vector dtscale{ getTimeStepScale(*fes.GetMesh()) };
    double rmin{ getJacobiGQ_RMin(fes.FEColl()->GetOrder()) };
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
        base_dt = calcMeshTimeStep(serialfes, 2.0, mfem::Element::Type::QUADRILATERAL) * opts.cfl;
    }
    else if (dimension == 3) {
        base_dt = calcMeshTimeStep(serialfes, 6.0, mfem::Element::Type::HEXAHEDRON) * opts.cfl;
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

void Solver::writeSimulationStatistics(const Time runtime){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::filesystem::path simExpPath("Exports/" + getRunModeTag() + "/" + model_.meshName_ + "/SimulationStats/");

    std::filesystem::create_directories(simExpPath);

    std::string path = (simExpPath / ("statistics_rank" + std::to_string(rank) + ".dat")).string();

    std::ofstream myfile(path, std::ios::app);
    if (myfile.is_open()) {
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
        myfile << "Memory Consumption (B): " << getCurrentMemoryUsage() << "\n";
        myfile.close();
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

void Solver::run()
{

    auto runStartTime = std::chrono::steady_clock::now();

#ifdef SHOW_TIMER_INFORMATION
    auto lastPrintTime{ std::chrono::steady_clock::now() };
    if (Mpi::WorldRank() == 0){
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "-------------SOLVER RUN INFORMATION-------------" << std::endl;
        printSimulationInformation(time_, dt_, opts_.final_time);
    }
#endif

    while (time_ <= opts_.final_time - 1e-8*dt_) {
        step();

#ifdef SHOW_TIMER_INFORMATION
        auto currentTime = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>
            (currentTime - lastPrintTime).count() >= 30.0)
        {
            if (Mpi::WorldRank() == 0){
                if (this->fields_.getNorml2() > 1e20){
                    std::cout << "------------------------------------------------------------------------" << std::endl;
                    std::cout << "Simulation is potentially unstable, verify manually and lower time step." << std::endl;
                    std::cout << "------------------------------------------------------------------------" << std::endl;
                }
                printSimulationInformation(time_, dt_, opts_.final_time);
                lastPrintTime = currentTime;
            }
        }
#endif
    }

    writeSimulationStatistics(std::chrono::duration<double>(std::chrono::steady_clock::now() - runStartTime).count());

}

void Solver::step()
{
    double truedt{ std::min(dt_, opts_.final_time - time_) };

    if (globalEvol_cache_) {
        std::array<mfem::ParGridFunction, 3> e, h;
        for(int d = 0; d < 3; ++d) {
            e[d].SetSpace(fields_.get(E, d).FESpace());
            h[d].SetSpace(fields_.get(H, d).FESpace());
            e[d] = fields_.get(E, d);
            h[d] = fields_.get(H, d);
        }
        globalEvol_cache_->advanceSGBCs(time_, truedt, e, h);
    }

    odeSolver_->Step(fields_.allDOFs(), time_, truedt);
    probesManager_.updateProbes(time_);
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