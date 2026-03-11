#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

#include <memory>

namespace maxwell
{

SGBCWrapper::~SGBCWrapper() = default;

const auto num_of_field_components = 2;
const auto num_of_max_dimensions = 3;
const auto num_of_field_blocks = num_of_field_components * num_of_max_dimensions;
const auto num_of_ghost_segments_per_field_comp = 2;

void SGBCWrapper::setAllSolverFields(const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields)
{
    // Warning: This sets the ACTIVE fields. In multi-state mode, this only affects currently loaded state.
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            sgbc_solver_fields_->get(f,d) = fields.get(f,d);
        }
    }
}

GeomTagToInteriorBoundary buildIntBdrInfo(const SGBCBoundaries& bdrInfo)
{
    GeomTagToInteriorBoundary res;
    if (bdrInfo.first.isOn){
        res[3] = bdrInfo.first.bdrCond;
    }
    if (bdrInfo.second.isOn){
        res[4] = bdrInfo.second.bdrCond;
    }
    return res;
}

GeomTagToBoundary buildBdrInfo()
{
    GeomTagToBoundary res;
    res[1] = BdrCond::SMA;
    res[2] = BdrCond::SMA;
    return res;
}

mfem::Mesh buildSGBCMesh(const SGBCProperties& sbcp)
{
    auto mesh = mfem::Mesh::MakeCartesian1D(sbcp.num_of_segments + 2, sbcp.material_width + 2 * sbcp.material_width / sbcp.num_of_segments);
    mesh.AddBdrPoint(1, 3);
    mesh.AddBdrPoint(mesh.GetNV() - 2, 4);
    mesh.SetAttribute(0, 2);
    mesh.SetAttribute(mesh.GetNE() - 1, 2);
    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4}); 
    mesh.FinalizeMesh();
    return mesh;
}

Model buildSGBCModel(mfem::Mesh& mesh, int* partitioning, const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo)
{
    Material vacuum = buildVacuumMaterial();
    GeomTagToMaterial geom_tag_sgbc_mat{{1, sbcp.material}, {2, vacuum}};
    GeomTagToInteriorBoundary gt2ib = buildIntBdrInfo(intBdrInfo);
    GeomTagToBoundary gt2b = buildBdrInfo();
    GeomTagToBoundaryInfo gtbdr(gt2b, gt2ib);
    return Model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), gtbdr, partitioning);
}

std::unique_ptr<SGBCWrapper> SGBCWrapper::buildSGBCWrapper(const SGBCProperties& sbcp)
{
    SGBCBoundaries bdrInfo = sbcp.sgbc_bdr_info;
    return std::unique_ptr<SGBCWrapper>(new SGBCWrapper(sbcp, bdrInfo));
}

std::unique_ptr<SGBCWrapper> SGBCWrapper::buildSGBCWrapperWithPEC(const SGBCProperties& sbcp)
{
    SGBCBoundaries bdrInfo;
    bdrInfo.first.bdrCond = BdrCond::PEC;
    bdrInfo.first.isOn = true;
    bdrInfo.second.bdrCond = BdrCond::PEC;
    bdrInfo.second.isOn = true;
    return std::unique_ptr<SGBCWrapper>(new SGBCWrapper(sbcp, bdrInfo));
}

SolverOptions buildSGBCSolverOptions(const SGBCProperties& sbcp)
{
    SolverOptions res;
    res.setOrder(sbcp.order);
    res.setUpwindAlpha(1.0);
    res.setODEType(ode_type::SDIRK33); 
    return res;
}

// [ADDED] Accessor
int SGBCWrapper::getStateSize() const {
    return solver_->getFields().allDOFs().Size();
}

// [ADDED] Context Switching
void SGBCWrapper::loadState(const SGBCState& state) {
    // Overwrite solver fields with the saved state
    solver_->getFields().allDOFs() = state.fields_state;
}

void SGBCWrapper::saveState(SGBCState& state) {
    // Save solver fields back to state
    state.fields_state = solver_->getFields().allDOFs();
}

void SGBCWrapper::updateFieldsWithGlobal(const std::array<mfem::ParGridFunction, 3>& e,
                                         const std::array<mfem::ParGridFunction, 3>& h,
                                         const SGBCState& context)
{
    const auto& dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    const NodePair& pair = context.global_pair;
    const int order_p1 = sbcp_.order + 1;
    const bool has_right = (pair.second != -1);
    const int right_dof_offset = dof_per_field_comp - 1;

    for (auto d : {X, Y, Z}){
        for (auto dof = 0; dof < order_p1; dof++){
            // Direct copy with bounds checking for safety-critical field access
            solver_->setFieldValue(FieldType::E, d, dof, e.at(d)[pair.first]);
            solver_->setFieldValue(FieldType::H, d, dof, h.at(d)[pair.first]);
            if (has_right) {
                solver_->setFieldValue(FieldType::E, d, right_dof_offset - dof, e.at(d)[pair.second]);
                solver_->setFieldValue(FieldType::H, d, right_dof_offset - dof, h.at(d)[pair.second]);
            }
        }
    }
}

void SGBCWrapper::getSGBCFields(const Array<int>& sub_to_global, const SGBCState& context, FieldGridFuncs& out)
{
    const auto left_ghost_border_dof = this->getProperties().order + 1;
    const auto local_field_size = this->solver_->getConstField(FieldType::E, X).Size();

    const auto first_idx = sub_to_global.Find(context.global_pair.first);
    int second_idx = -1;
    if (context.global_pair.second != -1){
        second_idx = sub_to_global.Find(context.global_pair.second);
    }

    if (first_idx == -1) return;

    int idx_left = left_ghost_border_dof;
    int idx_right = (local_field_size - 1) - left_ghost_border_dof;

    const int E_base = 0 * local_field_size;      // E fields at indices 0-2
    const int H_base = 3 * local_field_size;       // H fields at indices 3-5 (3 components each)

    for (auto f : {E, H}) {
        int base = (f == E) ? E_base : H_base;
        for (auto d : {X, Y, Z}) {
            int offset = base + d * local_field_size;
            out[f][d][first_idx] = context.fields_state[offset + idx_left];
            if (second_idx != -1) {
                out[f][d][second_idx] = context.fields_state[offset + idx_right];
            }
        }
    }
}

void SGBCWrapper::solve(const Time t, const Time dt)
{
    if (std::abs(dt) < 1e-9){
        return;
    }

    if (!temporal_warning_printed_) {
        temporal_warning_printed_ = true;

        double eps_r = sbcp_.material.getPermittivity();
        double mu_r = sbcp_.material.getPermeability();
        double sigma_solver = sbcp_.material.getConductivity();
        double dx = sbcp_.material_width / sbcp_.num_of_segments;

        // Calculate temporal parameters once
        // CRITICAL: Account for speed of light - convert spatial distance to temporal duration
        // crossing_time_SI [s] = dx [m] / v [m/s] where v = c / sqrt(eps_r*mu_r)
        // Therefore: crossing_time = dx * sqrt(eps_r*mu_r) / c_SI
        constexpr double c_si = physicalConstants::speedOfLight_SI;
        double crossing_time = (dx * std::sqrt(eps_r * mu_r)) / c_si;

        // Convert simulator time to SI seconds for proper comparison
        // Simulator uses normalized time: t_sim = t_si * c_si
        // Therefore: t_si = t_sim / c_si
        double dt_si = dt / c_si;

        // CRITICAL FIX: Proper dimensional relaxation time calculation
        // Physical relaxation time: τ = ε₀εᵣ / σ_SI [seconds]
        // In solver units where Z₀ = √(μ₀/ε₀), conductivity is stored as σ_solver = σ_SI × Z₀
        // Therefore to recover proper τ we must divide by Z₀ when computing in normalized units
        constexpr double z0_si = physicalConstants::freeSpaceImpedance_SI;
        constexpr double eps0_si = physicalConstants::vacuumPermittivity_SI;

        double relaxation_time = std::numeric_limits<double>::max();
        double loss_tangent = 0.0;  // tan(δ) = σ/(ωε) - dimensionless ratio

        if (sigma_solver > 0.0) {
            // Corrected formula for relaxation time in normalized units
            // tau = eps_r * eps0 / (sigma_solver / Z0) = eps_r * eps0 * Z0 / sigma_solver
            double sigma_si_recovered = sigma_solver / z0_si;
            relaxation_time = (eps_r * eps0_si) / sigma_si_recovered;

            // Material property characterization: loss tangent for adaptive factors
            // Estimate angular frequency: ω ~ 2π × max_freq (from solve context)
            // Loss tangent: tan(δ) = σ/(ω·ε) characterizes material behavior
            // This value is used to select appropriate time-stepping factors below
            loss_tangent = sigma_si_recovered / (eps_r * eps0_si);  // Normalized estimate
        }

        // Material-dependent temporal resolution factor (physics-based)
        // Metals (loss_tangent >> 1): Rapid charge relaxation, can use aggressive /10-15
        // Semiconductors (loss_tangent ~ 1): Intermediate case, /20-30
        // Weakly conducting (loss_tangent << 1): Slow decay, need conservative /50-100
        // Non-conducting (loss_tangent ≈ 0): No conductivity, skip relaxation criterion
        double relaxation_factor = 50.0;  // Default for weakly conducting
        if (loss_tangent > 10.0) {
            relaxation_factor = 10.0;     // Metals: fast relaxation, can be aggressive
        } else if (loss_tangent > 1.0) {
            relaxation_factor = 20.0;     // Semiconductors: intermediate
        } else if (loss_tangent > 0.1) {
            relaxation_factor = 50.0;     // Weakly conducting: conservative (default)
        }
        // else: non-conducting, relaxation_time already set to inf, factor not used

        double recommended_dt = crossing_time * 0.5;  // Wave criterion
        if (sigma_solver > 0.0) {
            recommended_dt = std::min(recommended_dt, relaxation_time / relaxation_factor);
        }

        bool is_coarse = (dt_si > recommended_dt);

        if (Mpi::WorldRank() == 0) {
            std::cout << "\n========================================================" << std::endl;
            if (is_coarse) {
                std::cout << "[SEVERE WARNING] SGBC Temporal Resolution is too coarse!" << std::endl;
            } else {
                std::cout << "[OK] SGBC Temporal Resolution is within good parameters!" << std::endl;
            }
            std::cout << "========================================================" << std::endl;
            std::cout << "  Current dt (SI)      : " << dt_si << " seconds" << std::endl;
            std::cout << "  Segment Crossing     : " << crossing_time << " seconds" << std::endl;
            if (sigma_solver > 0.0) {
                std::cout << "  Relaxation Time      : " << relaxation_time << " seconds" << std::endl;
                std::cout << "  Loss Tangent (tan δ) : " << loss_tangent << std::endl;
                std::cout << "  Relaxation Factor    : 1/" << relaxation_factor << std::endl;
            }

            if (is_coarse) {
                std::cout << "  Recommended dt (SI)  : < " << recommended_dt << " seconds" << std::endl;
                std::cout << "\n  The relaxation/crossing gradient is going to get mathematically smeared!" << std::endl;
                std::cout << "  Your implicit solver will remain stable, but your reflection physics" << std::endl;
                std::cout << "  will be wrong. You should lower your time step!" << std::endl;
            } else {
                std::cout << "  Max safe dt (SI)     : " << recommended_dt << " seconds" << std::endl;
                std::cout << "  Safety margin        : " << (recommended_dt / dt_si) << "x" << std::endl;
                std::cout << "\n  SGBC solver respects relaxation time material dynamics." << std::endl;
            }
            std::cout << "========================================================\n" << std::endl;
        }
    }

    this->solver_->setTime(t);
    this->solver_->getEvolTDO()->SetTime(t);

    this->solver_->getSolverOptions().setFinalTime(t + dt);
    this->solver_->getSolverOptions().setTimeStep(dt);

    this->solver_->setTimeStep(dt);

    this->solver_->step();
}

void checkSkinDepthResolution(const SGBCProperties& props)
{
    auto sigma_solver = props.material.getConductivity();

    if (sigma_solver <= 0.0){
        return;
    }

    constexpr auto Z0 = physicalConstants::freeSpaceImpedance_SI;
    constexpr auto mu0 = physicalConstants::vacuumPermeability_SI;
    constexpr auto eps0 = physicalConstants::vacuumPermittivity_SI;

    auto sigma_si = sigma_solver / Z0;
    auto mu_si = props.material.getPermeability() * mu0;
    auto eps_si = props.material.getPermittivity() * eps0;

    auto dx = props.material_width / static_cast<double>(props.num_of_segments);

    // Frequency at which skin depth equals element size
    auto f_max = 1.0 / (M_PI * mu_si * sigma_si * dx * dx);

    // Loss tangent for material characterization
    auto loss_tangent = sigma_si / (2.0 * M_PI * 1e9 * eps_si);  // At 1 GHz reference

    std::cout << "\n--- SGBC Spatial Resolution Check ---" << std::endl;
    std::cout << "Element dx                   : " << dx * 1e6 << " μm" << std::endl;
    std::cout << "Skin depth (at 1 GHz)        : " << (1.0 / sqrt(M_PI * 1e9 * mu_si * sigma_si) * 1e6) << " μm" << std::endl;
    std::cout << "Loss Tangent (tan δ @ 1 GHz) : " << loss_tangent << std::endl;
    std::cout << "Max resolved frequency (δ=dx): ";

    if (f_max > 1e10) {
        std::cout << f_max / 1e9 << " GHz" << std::endl;
        std::cout << "[EXCELLENT] Mesh resolution is outstanding for this material." << std::endl;
    }
    else if (f_max > 1e9) {
        std::cout << f_max / 1e9 << " GHz" << std::endl;
        std::cout << "[INFO] Mesh is well-resolved for standard RF/Microwave." << std::endl;
    }
    else if (f_max > 1e6) {
        std::cout << f_max / 1e6 << " MHz" << std::endl;
        std::cout << "[WARNING] Mesh may cause numerical tunneling for GHz pulses." << std::endl;
        std::cout << "          Consider increasing 'num_of_segments'." << std::endl;
    }
    else {
        std::cout << f_max / 1e3 << " kHz" << std::endl;
        std::cout << "[SEVERE] Mesh is severely coarse for this conductivity!" << std::endl;
        std::cout << "         Extreme numerical tunneling is highly likely." << std::endl;
    }
    std::cout << "----------------------------------\n" << std::endl;
}

SGBCWrapper::SGBCWrapper(const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo) :
sbcp_(sbcp)
{ 
    checkSkinDepthResolution(sbcp_);

    auto mesh = buildSGBCMesh(sbcp_);
    int* partitioning = mesh.GeneratePartitioning(1);
    
    Model model = buildSGBCModel(mesh, partitioning, sbcp_, intBdrInfo);
    Probes probes;
    // probes.exporterProbes.resize(1);
    // ExporterProbe ep;
    // ep.name = "InsideSGBC";
    // ep.visSteps = 10000;
    // probes.exporterProbes.at(0) = ep;
    Sources sources;
    SolverOptions opts = buildSGBCSolverOptions(sbcp_);
    opts.setExportEO(true);
    std::cout << "Assembling SGBC Solvers: " << std::endl;

    solver_ = std::make_unique<Solver>(model, probes, sources, opts);

    this->old_t_ = 0.0;
}

}