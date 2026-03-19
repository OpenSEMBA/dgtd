#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

#include <memory>
#include <algorithm>

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

mfem::Mesh buildSGBCMesh(const SGBCProperties& sbcp, int n_ghost)
{
    size_t total_segments = sbcp.totalSegments();
    double total_width = sbcp.totalWidth();
    double ghost_dx = total_width / total_segments;

    // Build vertex coordinates:
    //   n_ghost ghost elements | layer0 | ... | layerN | n_ghost ghost elements
    std::vector<double> vertices;

    // Left ghost vertices
    for (int i = 0; i <= n_ghost; ++i) {
        vertices.push_back(i * ghost_dx);
    }

    // Physical layer vertices (first vertex already placed above)
    double x = n_ghost * ghost_dx;
    for (const auto& layer : sbcp.layers) {
        double dx = layer.width / layer.num_of_segments;
        for (size_t i = 0; i < layer.num_of_segments; ++i) {
            x += dx;
            vertices.push_back(x);
        }
    }

    // Right ghost vertices
    for (int i = 1; i <= n_ghost; ++i) {
        vertices.push_back(x + i * ghost_dx);
    }

    int num_elements = static_cast<int>(vertices.size()) - 1;
    mfem::Mesh mesh(1, static_cast<int>(vertices.size()), num_elements, 0, 1);

    for (size_t i = 0; i < vertices.size(); ++i) {
        mesh.AddVertex(vertices[i]);
    }
    for (int i = 0; i < num_elements; ++i) {
        mesh.AddSegment(i, i + 1, 1);
    }
    mesh.AddBdrPoint(0, 1);
    mesh.AddBdrPoint(static_cast<int>(vertices.size()) - 1, 2);
    mesh.FinalizeMesh();

    // Interior boundary markers at ghost/material interfaces
    mesh.AddBdrPoint(n_ghost, 3);
    mesh.AddBdrPoint(n_ghost + static_cast<int>(total_segments), 4);

    // Assign element attributes:
    //   Elements 0..n_ghost-1:     left ghosts  -> vacuum_attr
    //   Physical layer elements:   attr = layer_index + 1
    //   Last n_ghost elements:     right ghosts -> vacuum_attr
    int vacuum_attr = static_cast<int>(sbcp.layers.size()) + 1;
    for (int i = 0; i < n_ghost; ++i) {
        mesh.SetAttribute(i, vacuum_attr);
    }

    int elem_idx = n_ghost;
    for (size_t li = 0; li < sbcp.layers.size(); ++li) {
        for (size_t s = 0; s < sbcp.layers[li].num_of_segments; ++s) {
            mesh.SetAttribute(elem_idx, static_cast<int>(li) + 1);
            elem_idx++;
        }
    }

    for (int i = 0; i < n_ghost; ++i) {
        mesh.SetAttribute(elem_idx + i, vacuum_attr);
    }

    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4});
    mesh.FinalizeMesh();
    return mesh;
}

Model buildSGBCModel(mfem::Mesh& mesh, int* partitioning, const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo)
{
    Material vacuum = buildVacuumMaterial();
    int vacuum_attr = static_cast<int>(sbcp.layers.size()) + 1;

    GeomTagToMaterial geom_tag_sgbc_mat;
    for (size_t li = 0; li < sbcp.layers.size(); ++li) {
        geom_tag_sgbc_mat.insert({static_cast<int>(li) + 1, sbcp.layers[li].material});
    }
    geom_tag_sgbc_mat.insert({vacuum_attr, vacuum});

    GeomTagToInteriorBoundary gt2ib = buildIntBdrInfo(intBdrInfo);
    GeomTagToBoundary gt2b = buildBdrInfo();
    GeomTagToBoundaryInfo gtbdr(gt2b, gt2ib);
    return Model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), gtbdr, partitioning, MPI_COMM_SELF);
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

std::unique_ptr<SGBCWrapper> SGBCWrapper::clone() const
{
    return buildSGBCWrapper(sbcp_);
}

SolverOptions buildSGBCSolverOptions(const SGBCProperties& sbcp)
{
    SolverOptions res;
    res.setOrder(sbcp.maxOrder());
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

int SGBCWrapper::getLocalFieldSize() const {
    return solver_->getConstField(FieldType::E, X).Size();
}

int SGBCWrapper::getLeftInterfaceIndex() const {
    return n_ghost_elements_ * (sbcp_.maxOrder() + 1);
}

int SGBCWrapper::getRightInterfaceIndex() const {
    int local_size = solver_->getConstField(FieldType::E, X).Size();
    return (local_size - 1) - n_ghost_elements_ * (sbcp_.maxOrder() + 1);
}

void SGBCWrapper::updateFieldsWithGlobal(const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields,
                                         const SGBCState& context)
{
    const auto& dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    const NodePair& pair = context.global_pair;
    const int total_ghost_dofs = n_ghost_elements_ * (sbcp_.maxOrder() + 1);
    const bool has_right = (pair.second != -1);
    const int right_dof_offset = dof_per_field_comp - 1;

    for (auto d : {X, Y, Z}){
        for (auto dof = 0; dof < total_ghost_dofs; dof++){
            solver_->setFieldValue(FieldType::E, d, dof, fields.get(E, d)[pair.first]);
            solver_->setFieldValue(FieldType::H, d, dof, fields.get(H, d)[pair.first]);
            if (has_right) {
                solver_->setFieldValue(FieldType::E, d, right_dof_offset - dof, fields.get(E, d)[pair.second]);
                solver_->setFieldValue(FieldType::H, d, right_dof_offset - dof, fields.get(H, d)[pair.second]);
            }
        }
    }
}

void SGBCWrapper::updateFieldsWithGlobalVector(const mfem::Vector& in, int ndofs, const SGBCState& context)
{
    const auto& dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    const NodePair& pair = context.global_pair;
    const int total_ghost_dofs = n_ghost_elements_ * (sbcp_.maxOrder() + 1);
    const bool has_right = (pair.second != -1);
    const int right_dof_offset = dof_per_field_comp - 1;

    for (int d = X; d <= Z; ++d){
        double e_left = in[d * ndofs + pair.first];
        double h_left = in[(3 + d) * ndofs + pair.first];
        double e_right = has_right ? in[d * ndofs + pair.second] : 0.0;
        double h_right = has_right ? in[(3 + d) * ndofs + pair.second] : 0.0;
        for (auto dof = 0; dof < total_ghost_dofs; dof++){
            solver_->setFieldValue(FieldType::E, d, dof, e_left);
            solver_->setFieldValue(FieldType::H, d, dof, h_left);
            if (has_right) {
                solver_->setFieldValue(FieldType::E, d, right_dof_offset - dof, e_right);
                solver_->setFieldValue(FieldType::H, d, right_dof_offset - dof, h_right);
            }
        }
    }
}

void SGBCWrapper::getSGBCFields(const Array<int>& sub_to_global, const SGBCState& context, SGBCHelperFields& out)
{
    const auto local_field_size = this->solver_->getConstField(FieldType::E, X).Size();

    const auto first_idx = sub_to_global.Find(context.global_pair.first);
    int second_idx = -1;
    if (context.global_pair.second != -1){
        second_idx = sub_to_global.Find(context.global_pair.second);
    }

    if (first_idx == -1) return;

    int idx_left = this->getLeftInterfaceIndex();
    int idx_right = this->getRightInterfaceIndex();

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

void SGBCWrapper::fillGlobalSGBCVec(const SGBCState& context, mfem::Vector& vec, int blockSize)
{
    const auto local_field_size = this->solver_->getConstField(FieldType::E, X).Size();

    const int idx_left  = this->getLeftInterfaceIndex();
    const int idx_right = this->getRightInterfaceIndex();

    const int gl = context.global_pair.first;
    const int gr = context.global_pair.second;

    for (auto f : {E, H}) {
        const int state_base = (f == E) ? 0 : 3 * local_field_size;
        for (auto d : {X, Y, Z}) {
            const int state_offset = state_base + d * local_field_size;
            vec[(f * 3 + d) * blockSize + gl] = context.fields_state[state_offset + idx_left];
            if (gr != -1) {
                vec[(f * 3 + d) * blockSize + gr] = context.fields_state[state_offset + idx_right];
            }
        }
    }
}

void SGBCWrapper::solve(const Time t, const Time dt)
{
    if (std::abs(dt) < 1e-9){
        return;
    }

    {
        static bool temporal_warning_printed = false;
        if (!temporal_warning_printed) {
            temporal_warning_printed = true;

            constexpr double c_si = physicalConstants::speedOfLight_SI;
            double dt_si = dt / c_si;
            double recommended_dt_si = recommended_dt_ / c_si;
            bool is_coarse = (dt_si > recommended_dt_si);

            if (Mpi::WorldRank() == 0) {
                std::cout << "\n========================================================" << std::endl;
                if (is_coarse) {
                    std::cout << "[SEVERE WARNING] SGBC Temporal Resolution is too coarse!" << std::endl;
                } else {
                    std::cout << "[OK] SGBC Temporal Resolution is within good parameters!" << std::endl;
                }
                std::cout << "========================================================" << std::endl;
                std::cout << "  Current dt (SI)      : " << dt_si << " seconds" << std::endl;
                std::cout << "  Recommended dt (SI)  : " << recommended_dt_si << " seconds" << std::endl;

                if (is_coarse) {
                    std::cout << "\n  Sub-stepping will be applied automatically." << std::endl;
                } else {
                    std::cout << "  Safety margin        : " << (recommended_dt_si / dt_si) << "x" << std::endl;
                }
                std::cout << "========================================================\n" << std::endl;
            }
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
    constexpr auto Z0 = physicalConstants::freeSpaceImpedance_SI;
    constexpr auto mu0 = physicalConstants::vacuumPermeability_SI;
    constexpr auto eps0 = physicalConstants::vacuumPermittivity_SI;

    for (size_t li = 0; li < props.layers.size(); ++li) {
        const auto& layer = props.layers[li];
        auto sigma_solver = layer.material.getConductivity();

        if (sigma_solver <= 0.0) continue;

        auto sigma_si = sigma_solver / Z0;
        auto mu_si = layer.material.getPermeability() * mu0;
        auto eps_si = layer.material.getPermittivity() * eps0;

        auto dx = layer.width / static_cast<double>(layer.num_of_segments);
        auto f_max = 1.0 / (M_PI * mu_si * sigma_si * dx * dx);
        auto loss_tangent = sigma_si / (2.0 * M_PI * 1e9 * eps_si);

        std::cout << "\n--- SGBC Spatial Resolution Check (Layer " << li + 1 << ") ---" << std::endl;
        std::cout << "Element dx                   : " << dx * 1e6 << " μm" << std::endl;
        std::cout << "Skin depth (at 1 GHz)        : " << (1.0 / sqrt(M_PI * 1e9 * mu_si * sigma_si) * 1e6) << " μm" << std::endl;
        std::cout << "Loss Tangent (tan δ @ 1 GHz) : " << loss_tangent << std::endl;
        std::cout << "Max resolved frequency (δ=dx): ";

        if (f_max > 1e10) {
            std::cout << f_max / 1e9 << " GHz" << std::endl;
            std::cout << "[EXCELLENT] Mesh resolution is outstanding for this material." << std::endl;
        } else if (f_max > 1e9) {
            std::cout << f_max / 1e9 << " GHz" << std::endl;
            std::cout << "[INFO] Mesh is well-resolved for standard RF/Microwave." << std::endl;
        } else if (f_max > 1e6) {
            std::cout << f_max / 1e6 << " MHz" << std::endl;
            std::cout << "[WARNING] Mesh may cause numerical tunneling for GHz pulses." << std::endl;
        } else {
            std::cout << f_max / 1e3 << " kHz" << std::endl;
            std::cout << "[SEVERE] Mesh is severely coarse for this conductivity!" << std::endl;
        }
        std::cout << "----------------------------------\n" << std::endl;
    }
}

SGBCWrapper::SGBCWrapper(const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo) :
sbcp_(sbcp),
n_ghost_elements_(std::max(3, static_cast<int>(sbcp.maxOrder()) + 1))
{ 
    checkSkinDepthResolution(sbcp_);

    auto mesh = buildSGBCMesh(sbcp_, n_ghost_elements_);
    int* partitioning = mesh.GeneratePartitioning(1);
    
    Model model = buildSGBCModel(mesh, partitioning, sbcp_, intBdrInfo);
    Probes probes;
    // probes.exporterProbes.resize(1);
    // ExporterProbe ep;
    // ep.name = "InsideSGBC";
    // ep.visSteps = 1000;
    // probes.exporterProbes.at(0) = ep;
    Sources sources;
    SolverOptions opts = buildSGBCSolverOptions(sbcp_);
    opts.setIsSGBCSolver(true);  // Mark as SGBC sub-solver to skip statistics
    std::cout << "Assembling SGBC Solvers: " << std::endl;

    solver_ = std::make_unique<Solver>(model, probes, sources, opts);

    this->old_t_ = 0.0;

    // Compute recommended_dt as the minimum across all layers.
    // CFL factor of 1.0: one wave crossing time per element.
    // This is the natural accuracy threshold for resolving wave propagation
    // through the layer. Works for all conductivities: low-sigma (hyperbolic,
    // need to resolve wavefronts) and high-sigma (parabolic, only skin-depth
    // attenuation matters — CFL 1.0 is conservative). The implicit solver
    // (BackwardEuler, L-stable) ensures unconditional stability regardless.
    // The dense direct solve makes sub-step cost negligible (O(n^2)
    // back-substitution on a small 1D system), so tighter CFL barely hurts.
    {
        constexpr double c_si = physicalConstants::speedOfLight_SI;

        recommended_dt_ = std::numeric_limits<double>::max();

        for (const auto& layer : sbcp_.layers) {
            double eps_r = layer.material.getPermittivity();
            double mu_r  = layer.material.getPermeability();
            double dx    = layer.width / layer.num_of_segments;

            double crossing_time = (dx * std::sqrt(eps_r * mu_r)) / c_si;
            double layer_dt = crossing_time * 0.5;

            recommended_dt_ = std::min(recommended_dt_, layer_dt);
        }

        // Convert from SI seconds to simulator time units
        recommended_dt_ *= c_si;
    }
}

}