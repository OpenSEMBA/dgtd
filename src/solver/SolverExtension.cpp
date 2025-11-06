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
    res[2] = (order + 1) * (num_of_segments + 2) - (order + 1); // end of mesh, SGBC element right boundary
    res[3] = (order + 1) * (num_of_segments + 2) - (order); // end of mesh, ghost element left boundary
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

void SBCSolver::update(const Time& dt)
{
    applyNodalToModalTransformation(); // S-1 x = q / Only flux rows belonging to the SGBC interfaces
    q_old_ = modal_values_;
    evol(q_old_, dt);
    applyModalToNodalTransformation(); // S q = x / Full matrix-vector operation
}

void SBCSolver::evol(const ModalValues& q_old, const Time& dt)
{
    for (auto i = 0; i < modal_values_.size(); i++){
        modal_values_[i] = std::exp(eigvals_[i] * dt) * q_old[i];
    }
}

void SBCSolver::applyNodalToModalTransformation()
{

    const auto field_offset = 3 * modal_values_.size() / 6;
    const auto dir_offset   =     modal_values_.size() / 6;

    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            modal_values_[f * field_offset + d * dir_offset + target_ids_[0]] = nodal_to_modal_rows_.at({f,d}).row_left_first.dot(nodal_values_);
            modal_values_[f * field_offset + d * dir_offset + target_ids_[1]] = nodal_to_modal_rows_.at({f,d}).row_left_second.dot(nodal_values_);
            modal_values_[f * field_offset + d * dir_offset + target_ids_[2]] = nodal_to_modal_rows_.at({f,d}).row_right_first.dot(nodal_values_);
            modal_values_[f * field_offset + d * dir_offset + target_ids_[3]] = nodal_to_modal_rows_.at({f,d}).row_right_second.dot(nodal_values_);
        }
    }
}

void SBCSolver::applyModalToNodalTransformation()
{
    auto temp = modal_to_nodal_matrix_ * modal_values_;
    if (temp.imag().cwiseAbs().sum() > 1e-5){
        throw std::runtime_error("Imaginary part over tolerance in ModalToNodalTransformation.");
    }
    nodal_values_ = temp.real();
}

void SBCSolver::initNodeIds(const std::pair<GlobalId, GlobalId>& el1_el2_ids)
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

SBCSolver::SBCSolver(const SBCProperties* sbcp, const std::pair<GlobalId, GlobalId>& global_dofs) : 
sbcp_(sbcp)
{
    const size_t number_of_field_components = 2;
    const size_t number_of_max_dimensions = 3;
    
    target_ids_ = buildTargetNodeIds(sbcp->order, sbcp->num_of_segments);
    initNodeIds(global_dofs);
    
    auto mesh  = Mesh::MakeCartesian1D(sbcp->num_of_segments + 2, sbcp->material_width + 2 * (sbcp->material_width / double(sbcp->num_of_segments)));
    auto pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    auto fec   = DG_FECollection(sbcp->order, 1, BasisType::GaussLobatto);
    auto pfes  = ParFiniteElementSpace(&pmesh, &fec);

    GeomTagToMaterial geom_tag_sgbc_mat{{1,sbcp_->material}};
    Model model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), GeomTagToBoundaryInfo{});
    auto global_operator = assembleGlobalOperator(model, pfes, sbcp->order);
    auto es = applyEigenSolverOnGlobalOperator(*global_operator);
    modal_to_nodal_matrix_ = es.eigenvectors();
    eigvals_ = es.eigenvalues(); // D = S-1 global_operator S, flattened to eigenvalues vector.

    modal_values_.resize(number_of_field_components * number_of_max_dimensions * target_ids_.size());
    modal_values_.setZero();

}

}