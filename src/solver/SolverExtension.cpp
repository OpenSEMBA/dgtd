#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

namespace maxwell
{

using namespace mfem;

SGBCWrapper::~SGBCWrapper() = default;

const auto num_of_field_components = 2;
const auto num_of_max_dimensions = 3;
const auto num_of_field_blocks = num_of_field_components * num_of_max_dimensions;
const auto num_of_ghost_segments_per_field_comp = 2;

std::vector<NodeId> buildTargetNodeIds(const size_t order, const size_t num_of_segments)
{
    std::vector<NodeId> res(4);
    res[0] = order; // start of mesh, ghost element right boundary   // |---Ghost---X-|-----SGBC---|--
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

Mesh buildSGBCMesh(const SGBCProperties& sbcp)
{
    auto mesh = Mesh::MakeCartesian1D(sbcp.num_of_segments + 2, sbcp.material_width + 2 * sbcp.material_width / sbcp.num_of_segments);
    mesh.AddBdrPoint(1, 3);
    mesh.AddBdrPoint(mesh.GetNV() - 2, 4);
    mesh.SetAttribute(0, 2);
    mesh.SetAttribute(mesh.GetNE() - 1, 2);
    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4}); // 1, 2 reserved for pure boundaries, 3, 4 reserved for interior boundaries.
    mesh.FinalizeMesh();
    return mesh;
}

Model buildSGBCModel(Mesh& mesh, int* partitioning, const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo)
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

SGBCWrapper::SGBCWrapper(const SGBCProperties& sbcp, const SGBCBoundaries& intBdrInfo) :
sbcp_(sbcp)
{ 
    
    auto mesh = buildSGBCMesh(sbcp_);
    int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
    
    Model model = buildSGBCModel(mesh, partitioning, sbcp_, intBdrInfo);
    Probes probes;
    Sources sources;
    SolverOptions opts;
    opts.setOrder(sbcp_.order);
    opts.setUpwindAlpha(1.0);
    opts.setODEType(ode_type::Trapezoidal); // Crank-Nicholson
    std::cout << "Assembling SGBC Solvers: " << std::endl;
    solver_ = std::make_unique<Solver>(model, probes, sources, opts);

}

}