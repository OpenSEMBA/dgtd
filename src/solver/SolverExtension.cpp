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
    res.setODEType(ode_type::SDIRK34); 
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

    for (auto d : {X, Y, Z}){ 
        for (auto dof = 0; dof < sbcp_.order + 1; dof++){ 
            
            // Direct Copy Logic (Stable)
            solver_->setFieldValue(FieldType::E, d, dof, e.at(d)[pair.first]); 
            solver_->setFieldValue(FieldType::H, d, dof, h.at(d)[pair.first]);
            if (context.global_pair.second != -1) { 
                solver_->setFieldValue(FieldType::E, d, dof_per_field_comp - 1 - dof, e.at(d)[pair.second]); 
                solver_->setFieldValue(FieldType::H, d, dof_per_field_comp - 1 - dof, h.at(d)[pair.second]);
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

    for (auto f : {E, H}) {
        for (auto d : {X, Y, Z}) {
            int block_idx = 0;
            if (f == E) {
                block_idx = d;
            } 
            else {
                block_idx = 3 + d;
            }   
            
            int offset = block_idx * local_field_size;

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

    this->solver_->setTime(t); 
    this->solver_->getEvolTDO()->SetTime(t);

    this->solver_->getSolverOptions().setFinalTime(t + dt);
    this->solver_->getSolverOptions().setTimeStep(dt); 
    
    this->solver_->setTimeStep(dt);

    this->solver_->step();
}

SGBCWrapper::SGBCWrapper(const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo) :
sbcp_(sbcp)
{ 
    auto mesh = buildSGBCMesh(sbcp_);
    int* partitioning = mesh.GeneratePartitioning(1);
    
    Model model = buildSGBCModel(mesh, partitioning, sbcp_, intBdrInfo);
    Probes probes;
    probes.exporterProbes.resize(1);
    ExporterProbe ep;
    ep.name = "InsideSGBC";
    ep.visSteps = 1000;
    probes.exporterProbes.at(0) = ep;
    Sources sources;
    SolverOptions opts = buildSGBCSolverOptions(sbcp_);
    opts.setExportEO(true);
    std::cout << "Assembling SGBC Solvers: " << std::endl;

    solver_ = std::make_unique<Solver>(model, probes, sources, opts);

    this->old_t_ = 0.0;
}

}