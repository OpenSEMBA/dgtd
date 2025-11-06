#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

namespace maxwell
{

using namespace mfem;

std::unique_ptr<SparseMatrix> assembleGlobalOperator(Model& model, ParFiniteElementSpace& pfes, size_t order)
{
    Probes pr;
    Sources src;
    EvolutionOptions eopts;
    eopts.order = order;
    ProblemDescription pd(model, pr, src, eopts);
    DGOperatorFactory<ParFiniteElementSpace> dgops(pd, pfes);
    return dgops.buildGlobalOperator();
}

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments)
{
    std::vector<NodeId> res(4);
    res[0] = order; // start of mesh, ghost element right boundary
    res[1] = order + 1; // start of mesh, SGBC element left boundary
    res[2] = (order + 1) * (num_of_segments + 2) - (order + 1) - 1; // end of mesh, SGBC element right boundary
    res[3] = (order + 1) * (num_of_segments + 2) - (order) - 1; // end of mesh, ghost element left boundary
    return res;
}

Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnGlobalOperator(const SparseMatrix& mat)
{
    MFEM_ASSERT(mat.Height() == mat.Width(), "Matrix must be square for EigenSolver.");

    using ESparseRow = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

    Eigen::Map<const ESparseRow> view(
        mat.Height(), mat.Width(), mat.NumNonZeroElems(),
        mat.GetI(), mat.GetJ(), mat.GetData());

    Eigen::MatrixXd M(view);

    Eigen::EigenSolver<Eigen::MatrixXd> res(M, true);
    MFEM_ASSERT(res.info() == Eigen::Success, "EigenSolver failed.");

    return res;
}

void SGBCSolver::setSGBCFieldValues(const SGBCNodalFields& in)
{
    int nodal_dofs = (sbcp_->num_of_segments + 2) * (sbcp_->order + 1);
    NodalValues nodal(6 * nodal_dofs);
    applyModalToNodalTransformation(nodal); // We perform S q = x of the previously stored q to obtain the latest updated x, if done after SGBCSolver initialization, this is simply zeroes.

    const auto field_offset = 3 * nodal_dofs;
    const auto dir_offset   =     nodal_dofs;

    for (auto f : {E, H}){ //Apply on ghost interface nodes.
        for (auto d : {X, Y, Z}){
            nodal[f * field_offset + d * dir_offset + dof_pair_.load.l_el1] = in.at(f).at(d).first;
            nodal[f * field_offset + d * dir_offset + dof_pair_.load.l_el2] = in.at(f).at(d).second;
        }
    }

    applyNodalToModalTransformation(nodal); // And now, with the updated nodal vector, we perform S-1 x = q and that q stays stored.

    // Technically we should be doing the thing at the bottom, extending to ghost dofs? For now this stays commented.

    // const auto& extra_dof = dof_pair_.load.l_el1; 

    // for (auto f : {E, H}){ //Extend to the rest of ghost element 1.
    //     for (auto d : {X, Y, Z}){
    //         for (auto i = 0; i < extra_dof; i++){
    //             nodal[f * field_offset + d * dir_offset + i] = vals.first;
    //         }
    //     }
    // }

    // for (auto f : {E, H}){ //Extend to the rest of ghost element 2.
    //     for (auto d : {X, Y, Z}){
    //         for (auto i = 0; i < extra_dof; i++){
    //             nodal[f * field_offset + d * dir_offset + dof_pair_.load.l_el2 + i] = vals.second;
    //         }
    //     }
    // }

}

SGBCNodalFields SGBCSolver::getSGBCFieldValues() const
{
    SGBCNodalFields res;

    int nodal_dofs = (sbcp_->num_of_segments + 2) * (sbcp_->order + 1);

    const auto field_offset = 3 * nodal_dofs;
    const auto dir_offset   =     nodal_dofs;

    NodalValues nodal(6 * nodal_dofs);
    applyModalToNodalTransformation(nodal); 

    for (auto f : {E, H}){ //Apply on ghost interface nodes.
        for (auto d : {X, Y, Z}){
            res.at(f).at(d).first = nodal[f * field_offset + d * dir_offset + dof_pair_.load.l_el1];
            res.at(f).at(d).second = nodal[f * field_offset + d * dir_offset + dof_pair_.load.l_el2];
        }
    }
    return res;
}

void SGBCSolver::update(const Time& dt)
{
    NodalValues temp;
    applyNodalToModalTransformation(temp); // S-1 x = q /
    q_old_ = modal_values_;
    evol(q_old_, dt);
    applyModalToNodalTransformation(temp); // S q = x /
}

void SGBCSolver::evol(const ModalValues& q_old, const Time& dt)
{
    for (auto i = 0; i < modal_values_.size(); i++){
        modal_values_[i] = std::exp(eigvals_[i] * dt) * q_old[i];
    }
}

void SGBCSolver::applyNodalToModalTransformation(const NodalValues& in)
{
    modal_values_ = nodal_to_modal_matrix_ * in;
}

void SGBCSolver::applyModalToNodalTransformation(NodalValues& out) const
{
    const Eigen::VectorXcd temp = modal_to_nodal_matrix_ * modal_values_;
    if (temp.imag().cwiseAbs().sum() > 1e-5){
        throw std::runtime_error("Imaginary part over tolerance in ModalToNodalTransformation.");
    }
    out = temp.real();
}

void SGBCSolver::initNodeIds(const std::pair<GlobalId, GlobalId>& el1_el2_ids)
{
    dof_pair_.load.g_el1 = el1_el2_ids.first;
    dof_pair_.load.g_el2 = el1_el2_ids.second;
    dof_pair_.load.l_el1 = target_ids_.front();
    dof_pair_.load.l_el2 = target_ids_.back();

    dof_pair_.unload.g_el1 = el1_el2_ids.first;
    dof_pair_.unload.g_el2 = el1_el2_ids.second;
    dof_pair_.unload.l_el1 = target_ids_.at(1);
    dof_pair_.unload.l_el2 = target_ids_.at(2);
}

SGBCSolver::SGBCSolver(const SGBCProperties* sbcp, const std::pair<GlobalId, GlobalId>& global_dofs) : 
sbcp_(sbcp)
{
    
    target_ids_ = buildTargetNodeIds(sbcp->order, sbcp->num_of_segments);
    initNodeIds(global_dofs);
    
    auto mesh  = Mesh::MakeCartesian1D(sbcp->num_of_segments + 2, sbcp->material_width + 2 * (sbcp->material_width / double(sbcp->num_of_segments)));
    int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
    auto pmesh = ParMesh(MPI_COMM_WORLD, mesh, partitioning);
    auto fec   = DG_FECollection(sbcp->order, 1, BasisType::GaussLobatto);
    auto pfes  = ParFiniteElementSpace(&pmesh, &fec);

    GeomTagToMaterial geom_tag_sgbc_mat{{1,sbcp_->material}};
    GeomTagToBoundary pecBdr{ {1,BdrCond::PEC}, {2,BdrCond::PEC} }; // For now, PEC is enforced on the boundaries
    Model model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), GeomTagToBoundaryInfo(pecBdr, GeomTagToInteriorBoundary()), partitioning);
    auto global_operator = assembleGlobalOperator(model, pfes, sbcp->order);
    std::cout << "Applying Eigen Solver on Global Operator." << std::endl;
    auto es = applyEigenSolverOnGlobalOperator(*global_operator);
    std::cout << "Eigen Solver success." << std::endl;
    modal_to_nodal_matrix_ = es.eigenvectors(); // S
    nodal_to_modal_matrix_ = modal_to_nodal_matrix_.inverse(); // S-1
    eigvals_ = es.eigenvalues(); // D = S-1 global_operator S, flattened to eigenvalues vector.

    modal_values_.resize(6 * pfes.GetNDofs());
    modal_values_.setZero();


}

}