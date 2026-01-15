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
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            sgbc_solver_fields_->get(f,d) = fields.get(f,d);
        }
    }
}

std::vector<NodeId> buildTargetNodeIds(const size_t order, const size_t num_of_segments)
{
    std::vector<NodeId> res(4);
    res[0] = order; // start of mesh, ghost element right boundary    // |---Ghost---X-|-----SGBC---|--
    res[1] = order + 1; // start of mesh, SGBC element left boundary // |---Ghost-----|-X---SGBC---|--
    res[2] = (order + 1) * (num_of_segments + 2) - (order + 1) - 1; // end of mesh, SGBC element right boundary // --|---SGBC---X-|-----Ghost---|
    res[3] = (order + 1) * (num_of_segments + 2) - (order) - 1; // end of mesh, ghost element left boundary     // --|---SGBC-----|-X---Ghost---|
    return res;
}

void SGBCWrapper::initNodeIds(const std::vector<NodeId>& target_ids) // To be redone
{
    dof_pair_.load_el1 = target_ids.front();
    dof_pair_.load_el2 = target_ids.back();

    dof_pair_.unload_el1 = target_ids.at(1);
    dof_pair_.unload_el2 = target_ids.at(2);
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
    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4}); // 1, 2 reserved for pure boundaries, 3, 4 reserved for interior boundaries.
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
    SGBCBoundaries bdrInfo;
    bdrInfo.first.isOn = false;
    bdrInfo.second.isOn = false;
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
    res.setODEType(ode_type::Trapezoidal); // Crank-Nicolson    
    return res;
}

void SGBCWrapper::updateFieldsWithGlobal(const std::array<mfem::ParGridFunction, 3>& e, const std::array<mfem::ParGridFunction, 3>& h, const NodePair& pair)
{
    const auto& dof_per_field_comp = solver_->getConstField(FieldType::E, X).Size();
    for (auto d : {X, Y, Z}){ 
        for (auto dof = 0; dof < sbcp_.order + 1; dof++){ 
            solver_->setFieldValue(FieldType::E, d, dof, e.at(d)[pair.first]); // Loading values on left side ghost
            solver_->setFieldValue(FieldType::H, d, dof, h.at(d)[pair.first]);
            solver_->setFieldValue(FieldType::E, d, dof_per_field_comp - 1 - dof, e.at(d)[pair.second]); // Loading values on right side ghost
            solver_->setFieldValue(FieldType::H, d, dof_per_field_comp - 1 - dof, h.at(d)[pair.second]);
        }
    }
}

void SGBCWrapper::getSGBCFields(const Array<int>& sub_to_global, const NodePair& pair, FieldGridFuncs& out)
{
    const auto left_ghost_border_dof = this->getProperties().order + 1;
    const auto local_field_size = this->solver_->getConstField(FieldType::E, X).Size();
    const auto first_idx = sub_to_global.Find(pair.first);
    const auto second_idx =sub_to_global.Find(pair.second);
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            out[f][d][first_idx] = this->solver_->getConstField(f,d)[left_ghost_border_dof];
            out[f][d][second_idx] = this->solver_->getConstField(f,d)[(local_field_size - 1) - left_ghost_border_dof];
        }
    }
}

void SGBCWrapper::solve(const Time t, const Time dt)
{
    if (std::abs(dt) < 1e-9){
        return;
    }
    this->solver_->getSolverOptions().setFinalTime(t + dt);
    this->solver_->getSolverOptions().setTimeStep(dt);
    this->solver_->getEvolTDO()->SetTime(t);
    this->solver_->step();
}

SGBCWrapper::SGBCWrapper(const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo) :
sbcp_(sbcp)
{ 
    
    auto mesh = buildSGBCMesh(sbcp_);
    int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
    
    Model model = buildSGBCModel(mesh, partitioning, sbcp_, intBdrInfo);
    Probes probes;
    // probes.exporterProbes.resize(1);
    // ExporterProbe ep;
    // ep.name = "InsideSGBC";
    // ep.visSteps = 1000;
    // probes.exporterProbes.at(0) = ep;
    Sources sources;
    SolverOptions opts = buildSGBCSolverOptions(sbcp_);
    opts.setExportEO(true);
    std::cout << "Assembling SGBC Solvers: " << std::endl;

    solver_ = std::make_unique<Solver>(model, probes, sources, opts);

    this->old_t_ = 0.0;

}

}