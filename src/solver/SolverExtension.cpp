#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

#include <cstring>
#include <algorithm>
#include <cmath>
#include <limits>

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

namespace {

GeomTagToMaterialInfo buildVacuumMaterialInfo(const mfem::ParSubMesh& mesh)
{
    GeomTagToMaterialInfo res;
    for (int e = 0; e < mesh.GetNE(); ++e) {
        res.gt2m.emplace(mesh.GetAttribute(e), buildVacuumMaterial());
    }
    return res;
}

const PMLRegionProperties& getCommonPMLProperties(const GeomTagToPMLRegion& pml_regions)
{
    if (pml_regions.empty()) {
        throw std::runtime_error("PML properties are not available for an empty volumetric PML region set.");
    }

    const auto& props = pml_regions.begin()->second;
    for (const auto& [tag, other] : pml_regions) {
        if (!(props == other)) {
            throw std::runtime_error("All volumetric PML tags must currently share the same properties.");
        }
    }
    return props;
}

double computeSigmaMaxForThickness(const double thickness, const PMLRegionProperties& props)
{
    if (thickness <= 0.0) {
        return 0.0;
    }
    return -(static_cast<double>(props.grading_order) + 1.0)
        * std::log(props.target_reflection)
        / (2.0 * thickness);
}

double evaluatePMLSigmaAtPosition(
    const mfem::Vector& position,
    const PMLBoxGeometry& geometry,
    const PMLRegionProperties& props,
    const mfem::Vector& sigma_max_minus,
    const mfem::Vector& sigma_max_plus)
{
    double sigma = 0.0;
    const int dim = geometry.inner_min.Size();
    for (int d = 0; d < dim; ++d) {
        if (!geometry.active_axes[d]) {
            continue;
        }

        if (position[d] < geometry.inner_min[d] && geometry.thickness_minus[d] > 0.0) {
            const double distance = geometry.inner_min[d] - position[d];
            const double xi = std::min(1.0, std::max(0.0, distance / geometry.thickness_minus[d]));
            sigma += sigma_max_minus[d] * std::pow(xi, static_cast<double>(props.grading_order));
        } else if (position[d] > geometry.inner_max[d] && geometry.thickness_plus[d] > 0.0) {
            const double distance = position[d] - geometry.inner_max[d];
            const double xi = std::min(1.0, std::max(0.0, distance / geometry.thickness_plus[d]));
            sigma += sigma_max_plus[d] * std::pow(xi, static_cast<double>(props.grading_order));
        }
    }
    return sigma;
}

std::pair<int, double> decodeSignedVdof(const int signed_vdof)
{
    if (signed_vdof >= 0) {
        return std::make_pair(signed_vdof, 1.0);
    }
    return std::make_pair(-1 - signed_vdof, -1.0);
}

double computeMinNormalElementSize(const mfem::ParSubMesh& mesh, const PMLBoxGeometry& geometry)
{
    double min_size = std::numeric_limits<double>::infinity();
    const int dim = geometry.inner_min.Size();
    for (int e = 0; e < mesh.GetNE(); ++e) {
        mfem::Array<int> vertices;
        mesh.GetElementVertices(e, vertices);
        mfem::Vector elem_min(dim);
        mfem::Vector elem_max(dim);
        for (int d = 0; d < dim; ++d) {
            elem_min[d] = std::numeric_limits<double>::infinity();
            elem_max[d] = -std::numeric_limits<double>::infinity();
        }

        for (int i = 0; i < vertices.Size(); ++i) {
            const double* coords = mesh.GetVertex(vertices[i]);
            for (int d = 0; d < dim; ++d) {
                elem_min[d] = std::min(elem_min[d], coords[d]);
                elem_max[d] = std::max(elem_max[d], coords[d]);
            }
        }

        for (int d = 0; d < dim; ++d) {
            if (!geometry.active_axes[d]) {
                continue;
            }

            const bool overlaps_axis_shell =
                elem_min[d] < geometry.inner_min[d] || elem_max[d] > geometry.inner_max[d];
            if (!overlaps_axis_shell) {
                continue;
            }

            const double extent = elem_max[d] - elem_min[d];
            if (extent > 0.0) {
                min_size = std::min(min_size, extent);
            }
        }
    }

    if (!std::isfinite(min_size)) {
        return 0.0;
    }
    return min_size;
}

} // namespace

PMLWrapper::PMLWrapper(Model& model, const mfem::ParFiniteElementSpace& parent_fes)
    : geometry_(model.inferPMLBoxGeometry()),
      submesher_(model.getConstMesh(), model.buildPMLVolumeMarker())
{
    auto* submesh = submesher_.getParSubMesh();
    if (submesh == nullptr) {
        throw std::runtime_error("PMLWrapper failed to build a volumetric PML ParSubMesh.");
    }

    sub_fes_ = std::make_unique<mfem::ParFiniteElementSpace>(submesh, parent_fes.FEColl());
    mfem::SubMeshUtils::BuildVdofToVdofMap(*sub_fes_,
                                          parent_fes,
                                          submesh->GetFrom(),
                                          submesh->GetParentElementIDMap(),
                                          sub_to_parent_vdofs_);

    const auto& props = getCommonPMLProperties(model.getPMLRegions());
    const int dim = geometry_.inner_min.Size();
    sigma_max_minus_.SetSize(dim);
    sigma_max_plus_.SetSize(dim);
    for (int d = 0; d < dim; ++d) {
        sigma_max_minus_[d] = computeSigmaMaxForThickness(geometry_.thickness_minus[d], props);
        sigma_max_plus_[d] = computeSigmaMaxForThickness(geometry_.thickness_plus[d], props);
        sigma_max_ = std::max(sigma_max_, std::max(sigma_max_minus_[d], sigma_max_plus_[d]));
    }

    auto sigma_function = [this, &props](const mfem::Vector& pos) {
        return evaluatePMLSigmaAtPosition(pos, geometry_, props, sigma_max_minus_, sigma_max_plus_);
    };
    mfem::FunctionCoefficient sigma_coeff(sigma_function);
    mfem::ParGridFunction sigma_gf(sub_fes_.get());
    sigma_gf.UseDevice(true);
    sigma_gf.ProjectCoefficient(sigma_coeff);
    sigma_profile_.SetSize(sub_fes_->GetNDofs());
    sigma_profile_.UseDevice(true);
    sigma_profile_ = sigma_gf;
    min_normal_element_size_ = computeMinNormalElementSize(*submesh, geometry_);

    Model pml_model(*submesh, buildVacuumMaterialInfo(*submesh));
    Probes probes;
    Sources sources;
    EvolutionOptions local_opts;
    ProblemDescription pd(pml_model, probes, sources, local_opts);
    DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, *sub_fes_);
    matched_layer_operator_ = dgops.buildMatchedConductiveOperator<mfem::ParBilinearForm>(sigma_coeff);
}

void PMLWrapper::gatherParentFields(
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& parent_fields,
    mfem::Vector& sub_state) const
{
    const int ndofs = sub_fes_->GetNDofs();
    sub_state.SetSize(getStateSize());
    sub_state.UseDevice(true);
    for (int comp = 0; comp < 6; ++comp) {
        const FieldType field = comp < 3 ? FieldType::E : FieldType::H;
        const Direction dir = static_cast<Direction>(comp % 3);
        for (int vdof = 0; vdof < ndofs; ++vdof) {
            auto [parent_vdof, sign] = decodeSignedVdof(sub_to_parent_vdofs_[vdof]);
            sub_state[comp * ndofs + vdof] = sign * parent_fields.get(field, dir)[parent_vdof];
        }
    }
}

void PMLWrapper::scatterStateToParent(
    const mfem::Vector& sub_state,
    Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& parent_fields) const
{
    const int ndofs = sub_fes_->GetNDofs();
    for (int comp = 0; comp < 6; ++comp) {
        const FieldType field = comp < 3 ? FieldType::E : FieldType::H;
        const Direction dir = static_cast<Direction>(comp % 3);
        for (int vdof = 0; vdof < ndofs; ++vdof) {
            auto [parent_vdof, sign] = decodeSignedVdof(sub_to_parent_vdofs_[vdof]);
            parent_fields.get(field, dir)[parent_vdof] = sign * sub_state[comp * ndofs + vdof];
        }
    }
}

void PMLWrapper::initializeStateFromParent(
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& parent_fields,
    PMLRegionState& state) const
{
    if (state.auxiliary_state.Size() != getStateSize()) {
        state.init(getStateSize());
    }
    gatherParentFields(parent_fields, state.auxiliary_state);
}

void PMLWrapper::applyImplicitCorrection(
    double dt,
    Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& parent_fields,
    PMLRegionState& state) const
{
    if (state.auxiliary_state.Size() != getStateSize()) {
        state.init(getStateSize());
    }

    mfem::Vector rhs;
    rhs.UseDevice(true);
    gatherParentFields(parent_fields, rhs);
    state.auxiliary_state = rhs;

    if (dt <= 0.0 || matched_layer_operator_ == nullptr || sigma_max_ <= 0.0) {
        scatterStateToParent(state.auxiliary_state, parent_fields);
        return;
    }

    class PMLJacobian final : public mfem::Operator {
    public:
        PMLJacobian(const mfem::SparseMatrix& op, double dt)
            : mfem::Operator(op.Height()), op_(op), dt_(dt), tmp_(op.Height())
        {
            tmp_.UseDevice(true);
        }

        void Mult(const mfem::Vector& v, mfem::Vector& Jv) const override
        {
            op_.Mult(v, tmp_);
            Jv = v;
            Jv.Add(-dt_, tmp_);
        }

    private:
        const mfem::SparseMatrix& op_;
        double dt_;
        mutable mfem::Vector tmp_;
    } jacobian(*matched_layer_operator_, dt);

    mfem::GMRESSolver gmres;
    gmres.SetOperator(jacobian);
    gmres.SetRelTol(1e-8);
    gmres.SetMaxIter(200);
    gmres.SetKDim(50);
    gmres.SetPrintLevel(0);
    gmres.Mult(rhs, state.auxiliary_state);

    scatterStateToParent(state.auxiliary_state, parent_fields);
}

std::unique_ptr<SGBCWrapper> SGBCWrapper::buildSGBCWrapper(const SGBCProperties& sbcp, double simulation_final_time, const ExporterProbe* exporter_probe)
{
    SGBCBoundaries bdrInfo = sbcp.sgbc_bdr_info;
    return std::unique_ptr<SGBCWrapper>(new SGBCWrapper(sbcp, bdrInfo, simulation_final_time, exporter_probe));
}

std::unique_ptr<SGBCWrapper> SGBCWrapper::buildSGBCWrapperWithPEC(const SGBCProperties& sbcp, double simulation_final_time, const ExporterProbe* exporter_probe)
{
    SGBCBoundaries bdrInfo;
    bdrInfo.first.bdrCond = BdrCond::PEC;
    bdrInfo.first.isOn = true;
    bdrInfo.second.bdrCond = BdrCond::PEC;
    bdrInfo.second.isOn = true;
    return std::unique_ptr<SGBCWrapper>(new SGBCWrapper(sbcp, bdrInfo, simulation_final_time, exporter_probe));
}

std::unique_ptr<SGBCWrapper> SGBCWrapper::clone() const
{
    const ExporterProbe* exporter_probe = has_exporter_probe_ ? &exporter_probe_ : nullptr;
    return std::unique_ptr<SGBCWrapper>(new SGBCWrapper(sbcp_, intBdrInfo_, simulation_final_time_, exporter_probe));
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
    auto& dst = solver_->getFields().allDOFs();
    std::memcpy(dst.GetData(), state.fields_state.GetData(),
                dst.Size() * sizeof(double));
}

void SGBCWrapper::saveState(SGBCState& state) {
    const auto& src = solver_->getFields().allDOFs();
    std::memcpy(state.fields_state.GetData(), src.GetData(),
                src.Size() * sizeof(double));
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
    const int dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    const NodePair& pair = context.global_pair;
    const int total_ghost_dofs = n_ghost_elements_ * (sbcp_.maxOrder() + 1);
    const bool has_right = (pair.second != -1);
    const int right_dof_offset = dof_per_field_comp - 1;
    double* all = solver_->getFields().allDOFs().GetData();
    const double* R = context.rot;

    // Read global Cartesian fields at left DOF
    double eg[3], hg[3], el[3], hl[3];
    for (int d = 0; d < 3; ++d) {
        eg[d] = fields.get(E, static_cast<Direction>(d))[pair.first];
        hg[d] = fields.get(H, static_cast<Direction>(d))[pair.first];
    }
    // Rotate global → face-local: E_local = R * E_global
    for (int i = 0; i < 3; ++i) {
        el[i] = R[3*i]*eg[0] + R[3*i+1]*eg[1] + R[3*i+2]*eg[2];
        hl[i] = R[3*i]*hg[0] + R[3*i+1]*hg[1] + R[3*i+2]*hg[2];
    }
    for (int d = 0; d < 3; ++d) {
        double* e_block = all + d * dof_per_field_comp;
        double* h_block = all + (3 + d) * dof_per_field_comp;
        for (int dof = 0; dof < total_ghost_dofs; ++dof) {
            e_block[dof] = el[d];
            h_block[dof] = hl[d];
        }
    }
    if (has_right) {
        for (int d = 0; d < 3; ++d) {
            eg[d] = fields.get(E, static_cast<Direction>(d))[pair.second];
            hg[d] = fields.get(H, static_cast<Direction>(d))[pair.second];
        }
        for (int i = 0; i < 3; ++i) {
            el[i] = R[3*i]*eg[0] + R[3*i+1]*eg[1] + R[3*i+2]*eg[2];
            hl[i] = R[3*i]*hg[0] + R[3*i+1]*hg[1] + R[3*i+2]*hg[2];
        }
        for (int d = 0; d < 3; ++d) {
            double* e_block = all + d * dof_per_field_comp;
            double* h_block = all + (3 + d) * dof_per_field_comp;
            for (int dof = 0; dof < total_ghost_dofs; ++dof) {
                e_block[right_dof_offset - dof] = el[d];
                h_block[right_dof_offset - dof] = hl[d];
            }
        }
    }
}

void SGBCWrapper::updateFieldsWithGlobalVector(const mfem::Vector& in, int ndofs, const SGBCState& context)
{
    const int dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    const NodePair& pair = context.global_pair;
    const int total_ghost_dofs = n_ghost_elements_ * (sbcp_.maxOrder() + 1);
    const bool has_right = (pair.second != -1);
    const int right_dof_offset = dof_per_field_comp - 1;
    double* all = solver_->getFields().allDOFs().GetData();
    const double* R = context.rot;

    // Read global-frame fields at left DOF, rotate to face-local frame
    double eg[3], hg[3], el[3], hl[3];
    for (int d = 0; d < 3; ++d) {
        eg[d] = in[d * ndofs + pair.first];
        hg[d] = in[(3 + d) * ndofs + pair.first];
    }
    for (int i = 0; i < 3; ++i) {
        el[i] = R[3*i]*eg[0] + R[3*i+1]*eg[1] + R[3*i+2]*eg[2];
        hl[i] = R[3*i]*hg[0] + R[3*i+1]*hg[1] + R[3*i+2]*hg[2];
    }
    for (int d = 0; d < 3; ++d) {
        double* e_block = all + d * dof_per_field_comp;
        double* h_block = all + (3 + d) * dof_per_field_comp;
        for (int dof = 0; dof < total_ghost_dofs; ++dof) {
            e_block[dof] = el[d];
            h_block[dof] = hl[d];
        }
    }
    if (has_right) {
        for (int d = 0; d < 3; ++d) {
            eg[d] = in[d * ndofs + pair.second];
            hg[d] = in[(3 + d) * ndofs + pair.second];
        }
        for (int i = 0; i < 3; ++i) {
            el[i] = R[3*i]*eg[0] + R[3*i+1]*eg[1] + R[3*i+2]*eg[2];
            hl[i] = R[3*i]*hg[0] + R[3*i+1]*hg[1] + R[3*i+2]*hg[2];
        }
        for (int d = 0; d < 3; ++d) {
            double* e_block = all + d * dof_per_field_comp;
            double* h_block = all + (3 + d) * dof_per_field_comp;
            for (int dof = 0; dof < total_ghost_dofs; ++dof) {
                e_block[right_dof_offset - dof] = el[d];
                h_block[right_dof_offset - dof] = hl[d];
            }
        }
    }
}

void SGBCWrapper::updateFieldsWithInterpolatedGhost(double alpha, const SGBCState& context)
{
    const int dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    const int total_ghost_dofs = n_ghost_elements_ * (sbcp_.maxOrder() + 1);
    const bool has_right = (context.global_pair.second != -1);
    const int right_dof_offset = dof_per_field_comp - 1;
    double* all = solver_->getFields().allDOFs().GetData();

    const double beta = 1.0 - alpha;

    // Left ghost: interpolate face-local fields
    for (int d = 0; d < 3; ++d) {
        double ev = beta * context.ghost_old[d]     + alpha * context.ghost_new[d];
        double hv = beta * context.ghost_old[3 + d] + alpha * context.ghost_new[3 + d];
        double* e_block = all + d * dof_per_field_comp;
        double* h_block = all + (3 + d) * dof_per_field_comp;
        for (int dof = 0; dof < total_ghost_dofs; ++dof) {
            e_block[dof] = ev;
            h_block[dof] = hv;
        }
    }
    if (has_right) {
        for (int d = 0; d < 3; ++d) {
            double ev = beta * context.ghost_old[6 + d]     + alpha * context.ghost_new[6 + d];
            double hv = beta * context.ghost_old[6 + 3 + d] + alpha * context.ghost_new[6 + 3 + d];
            double* e_block = all + d * dof_per_field_comp;
            double* h_block = all + (3 + d) * dof_per_field_comp;
            for (int dof = 0; dof < total_ghost_dofs; ++dof) {
                e_block[right_dof_offset - dof] = ev;
                h_block[right_dof_offset - dof] = hv;
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
    const double* R = context.rot;

    // Left interface: read face-local fields and rotate to global
    double el[3], hl[3];
    for (int d = 0; d < 3; ++d) {
        el[d] = context.fields_state[d * local_field_size + idx_left];
        hl[d] = context.fields_state[(3 + d) * local_field_size + idx_left];
    }
    // R^T * local = global  (R is orthogonal)
    for (int d = 0; d < 3; ++d) {
        out[E][d][first_idx] = R[d]*el[0] + R[3+d]*el[1] + R[6+d]*el[2];
        out[H][d][first_idx] = R[d]*hl[0] + R[3+d]*hl[1] + R[6+d]*hl[2];
    }

    if (second_idx != -1) {
        for (int d = 0; d < 3; ++d) {
            el[d] = context.fields_state[d * local_field_size + idx_right];
            hl[d] = context.fields_state[(3 + d) * local_field_size + idx_right];
        }
        for (int d = 0; d < 3; ++d) {
            out[E][d][second_idx] = R[d]*el[0] + R[3+d]*el[1] + R[6+d]*el[2];
            out[H][d][second_idx] = R[d]*hl[0] + R[3+d]*hl[1] + R[6+d]*hl[2];
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
    const double* R = context.rot;

    // Left interface: read face-local fields and rotate to global
    double el[3], hl[3];
    for (int d = 0; d < 3; ++d) {
        el[d] = context.fields_state[d * local_field_size + idx_left];
        hl[d] = context.fields_state[(3 + d) * local_field_size + idx_left];
    }
    for (int d = 0; d < 3; ++d) {
        vec[d * blockSize + gl]       = R[d]*el[0] + R[3+d]*el[1] + R[6+d]*el[2];
        vec[(3 + d) * blockSize + gl] = R[d]*hl[0] + R[3+d]*hl[1] + R[6+d]*hl[2];
    }

    if (gr != -1) {
        for (int d = 0; d < 3; ++d) {
            el[d] = context.fields_state[d * local_field_size + idx_right];
            hl[d] = context.fields_state[(3 + d) * local_field_size + idx_right];
        }
        for (int d = 0; d < 3; ++d) {
            vec[d * blockSize + gr]       = R[d]*el[0] + R[3+d]*el[1] + R[6+d]*el[2];
            vec[(3 + d) * blockSize + gr] = R[d]*hl[0] + R[3+d]*hl[1] + R[6+d]*hl[2];
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
                if (is_coarse) {
                    std::cout << "[SGBC] WARNING: parent dt=" << dt_si*1e12 << " ps > recommended "
                              << recommended_dt_si*1e12 << " ps — sub-stepping applied.\n" << std::flush;
                } else {
                    std::cout << "[SGBC] Temporal OK: parent dt=" << dt_si*1e12 << " ps, recommended "
                              << recommended_dt_si*1e12 << " ps (margin " << (recommended_dt_si / dt_si) << "x)\n" << std::flush;
                }
            }
        }
    }

    this->solver_->setTime(t);
    this->solver_->getEvolTDO()->SetTime(t);

    this->solver_->setTimeStep(dt);

    this->solver_->step(false);
}

void SGBCWrapper::updateProbes(const Time t)
{
    this->solver_->setTime(t);
    this->solver_->getEvolTDO()->SetTime(t);
    this->solver_->updateProbes();
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

SGBCWrapper::SGBCWrapper(const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo, double simulation_final_time, const ExporterProbe* exporter_probe) :
sbcp_(sbcp),
intBdrInfo_(intBdrInfo),
simulation_final_time_(simulation_final_time),
has_exporter_probe_(exporter_probe != nullptr),
exporter_probe_(exporter_probe != nullptr ? *exporter_probe : ExporterProbe{}),
n_ghost_elements_(std::max(3, static_cast<int>(sbcp.maxOrder()) + 1))
{ 
    auto mesh = buildSGBCMesh(sbcp_, n_ghost_elements_);
    int* partitioning = mesh.GeneratePartitioning(1);
    
    Model model = buildSGBCModel(mesh, partitioning, sbcp_, intBdrInfo);
    Probes probes;
    if (has_exporter_probe_) {
        auto ep = exporter_probe_;
        const auto tag = sbcp_.geom_tags.empty() ? 0 : sbcp_.geom_tags.front();
        ep.name = "InsideSGBC_tag" + std::to_string(tag);
        probes.exporterProbes.push_back(ep);
    }
    Sources sources;
    SolverOptions opts = buildSGBCSolverOptions(sbcp_);
    opts.setFinalTime(simulation_final_time_);
    opts.setIsSGBCSolver(true);  // Mark as SGBC sub-solver to skip statistics

    solver_ = std::make_unique<Solver>(model, probes, sources, opts);

    this->old_t_ = 0.0;

    // Compute recommended_dt as the minimum across all layers.
    // Base CFL: half the wave crossing time per element.
    // For layers that are many skin depths thick (N_delta >> 1), the physics
    // is diffusion-dominated and the L-stable implicit solver can safely take
    // much larger steps. We relax the CFL proportionally to N_delta^2,
    // capped to avoid excessive jumps.
    {
        constexpr double c_si = physicalConstants::speedOfLight_SI;
        constexpr double cfl_relax_cap = 50.0;

        recommended_dt_ = std::numeric_limits<double>::max();

        for (const auto& layer : sbcp_.layers) {
            double eps_r = layer.material.getPermittivity();
            double mu_r  = layer.material.getPermeability();
            double dx    = layer.width / layer.num_of_segments;

            double crossing_time = (dx * std::sqrt(eps_r * mu_r)) / c_si;
            double layer_dt = crossing_time * 0.5;

            // Relax CFL for opaque layers: N_delta > 3 means wave is
            // heavily attenuated, so temporal resolution of wave transit
            // is unnecessary. The implicit solver handles the stiffness.
            double nd = layer.n_skin_depths;
            if (nd > 3.0) {
                double relax = std::min(nd * nd, cfl_relax_cap);
                layer_dt *= relax;
            }

            recommended_dt_ = std::min(recommended_dt_, layer_dt);
        }

        // Convert from SI seconds to simulator time units
        recommended_dt_ *= c_si;
    }

    // Print spatial + temporal resolution check once (first wrapper only).
    {
        static bool resolutionChecked = false;
        if (!resolutionChecked) {
            resolutionChecked = true;
            checkSkinDepthResolution(sbcp_);
            constexpr double c_si = physicalConstants::speedOfLight_SI;
            double rec_dt_si = recommended_dt_ / c_si;
            std::cout << "  SGBC recommended dt   : " << rec_dt_si * 1e12 << " ps"
                      << "  (CFL-relaxed for " << sbcp_.layers[0].n_skin_depths
                      << " skin depths)\n" << std::endl;
        }
    }
}

}