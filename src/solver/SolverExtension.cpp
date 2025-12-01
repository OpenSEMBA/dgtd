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
    res[0] = order; // start of mesh, ghost element right boundary   // |---Ghost---X-|-----SGBC---|--
    res[1] = order + 1; // start of mesh, SGBC element left boundary // |---Ghost-----|-X---SGBC---|--
    res[2] = (order + 1) * (num_of_segments + 2) - (order + 1) - 1; // end of mesh, SGBC element right boundary // --|---SGBC---X-|-----Ghost---|
    res[3] = (order + 1) * (num_of_segments + 2) - (order) - 1; // end of mesh, ghost element left boundary     // --|---SGBC-----|-X---Ghost---|
    return res;
}

Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnLocalOperator(const Eigen::MatrixXd& mat)
{
    MFEM_ASSERT(mat.cols() == mat.rows(), "Matrix must be square for EigenSolver.");

    Eigen::EigenSolver<Eigen::MatrixXd> res(mat, true);
    MFEM_ASSERT(res.info() == Eigen::Success, "EigenSolver failed.");

    return res;
}

// Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnLocalOperator(const SparseMatrix& mat)
// {
//     MFEM_ASSERT(mat.Height() == mat.Width(), "Matrix must be square for EigenSolver.");

//     using ESparseRow = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

//     Eigen::Map<const ESparseRow> view(
//         mat.Height(), mat.Width(), mat.NumNonZeroElems(),
//         mat.GetI(), mat.GetJ(), mat.GetData());

//     Eigen::MatrixXd M(view);

//     Eigen::EigenSolver<Eigen::MatrixXd> res(M, true);
//     MFEM_ASSERT(res.info() == Eigen::Success, "EigenSolver failed.");

//     return res;
// }

void SGBCSolver::setFullNodalState(const FullNodalFields& in)
{
    const auto& ndofs = in.at(0).at(0).Size();
    NodalValues nodal(6 * ndofs);
    
    const auto field_offset = 3 * ndofs;
    const auto dir_offset   =     ndofs;
    
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            for (auto v = 0; v < ndofs; v++){
                nodal[f * field_offset + d * dir_offset + v] = in.at(f).at(d)[v];
            }
        }
    }

    applyLocalNodalToLocalModalTransformation(nodal);

}

FullNodalFields SGBCSolver::getFullNodalState() const
{
    FullNodalFields res;

    int ndofs = (sbcp_->num_of_segments + 2) * (sbcp_->order + 1);
    const auto field_offset = 3 * ndofs;
    const auto dir_offset   =     ndofs;
    
    NodalValues nodal(6 * ndofs);
    applyModalToNodalTransformation(nodal); 
    
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            res.at(f).at(d).SetSize(ndofs);
            for (auto v = 0; v < ndofs; v++){
                res.at(f).at(d)[v] = nodal[f * field_offset + d * dir_offset + v];
            }
        }
    }

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
            nodal[f * field_offset + d * dir_offset + dof_pair_.load_el1] = in.at(f).at(d).first;
            nodal[f * field_offset + d * dir_offset + dof_pair_.load_el2] = in.at(f).at(d).second;
        }
    }

    applyLocalNodalToLocalModalTransformation(nodal); // And now, with the updated nodal vector, we perform S-1 x = q and that q stays stored.

    const auto& extra_dof = dof_pair_.load_el1; 

    for (auto f : {E, H}){ //Extend to the rest of ghost element 1.
        for (auto d : {X, Y, Z}){
            for (auto i = 0; i < extra_dof; i++){
                nodal[f * field_offset + d * dir_offset + i] = in.at(f).at(d).first;
            }
        }
    }

    for (auto f : {E, H}){ //Extend to the rest of ghost element 2.
        for (auto d : {X, Y, Z}){
            for (auto i = 0; i < extra_dof; i++){
                nodal[f * field_offset + d * dir_offset + dof_pair_.load_el2 + i] = in.at(f).at(d).second;
            }
        }
    }

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
            res.at(f).at(d).first = nodal[f * field_offset + d * dir_offset + dof_pair_.load_el1];
            res.at(f).at(d).second = nodal[f * field_offset + d * dir_offset + dof_pair_.load_el2];
        }
    }
    return res;
}

void SGBCSolver::update(const Time& dt)
{
    local_modal_values_old_ = local_modal_values_;
    evol(local_modal_values_old_, dt);
}

void SGBCSolver::evol(const ModalValues& q_old, const Time& dt)
{
    for (auto i = 0; i < local_modal_values_.size(); i++){
        local_modal_values_[i] = std::exp(lambda_[i] * dt) * q_old[i];
    }
}

void SGBCSolver::applyLocalNodalToLocalModalTransformation(const NodalValues& in)
{
    local_modal_values_ = Sinv_ * in;
}

void SGBCSolver::applyModalToNodalTransformation(NodalValues& out) const
{
    const Eigen::VectorXcd temp = S_ * local_modal_values_;
    if (temp.imag().cwiseAbs().sum() > 1e-5){
        throw std::runtime_error("Imaginary part over tolerance in ModalToNodalTransformation.");
    }
    out = temp.real();
}

void SGBCSolver::initNodeIds(const std::vector<NodeId>& target_ids) // To be redone
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

Mesh buildSGBCMesh(const SGBCProperties* sbcp)
{
    auto mesh = Mesh::MakeCartesian1D(sbcp->num_of_segments + 2, sbcp->material_width + 2 * sbcp->material_width / sbcp->num_of_segments);
    mesh.AddBdrPoint(1, 3);
    mesh.AddBdrPoint(mesh.GetNV() - 2, 4);
    mesh.SetAttribute(0, 2);
    mesh.SetAttribute(mesh.GetNE() - 1, 2);
    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4}); // 1, 2 reserved for pure boundaries, 3, 4 reserved for interior boundaries.
    mesh.FinalizeMesh();
    return mesh;
}

Model buildSGBCModel(Mesh& mesh, int* partitioning, const SGBCProperties* sbcp, const SGBCBoundaries& intBdrInfo)
{
    Material vacuum(1.0, 1.0, 0.0);
    GeomTagToMaterial geom_tag_sgbc_mat{{1, sbcp->material}, {2, vacuum}};
    GeomTagToInteriorBoundary gt2ib = buildIntBdrInfo(intBdrInfo);
    GeomTagToBoundary gt2b = buildBdrInfo();
    GeomTagToBoundaryInfo gtbdr(gt2b, gt2ib);
    return Model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), gtbdr, partitioning);
}

std::unique_ptr<SGBCSolver> SGBCSolver::buildSGBCSolver(const SGBCProperties* sbcp)
{
    SGBCBoundaries bdrInfo;
    bdrInfo.first.isOn = false;
    bdrInfo.second.isOn = false;
    return std::unique_ptr<SGBCSolver>(new SGBCSolver(sbcp, bdrInfo));
}

std::unique_ptr<SGBCSolver> SGBCSolver::buildSGBCSolverWithPEC(const SGBCProperties* sbcp)
{
    SGBCBoundaries bdrInfo;
    bdrInfo.first.bdrCond = BdrCond::PEC;
    bdrInfo.first.isOn = true;
    bdrInfo.second.bdrCond = BdrCond::PEC;
    bdrInfo.second.isOn = true;
    return std::unique_ptr<SGBCSolver>(new SGBCSolver(sbcp, bdrInfo));
}

size_t SGBCSolver::getLocalModalSize()
{
    return local_modal_values_.size();
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> buildLocalAndGhostMatrices(const SparseMatrix& global, const int order)
{
    const auto block_width = order + 1;
    const auto num_of_ghost_segments = 2;
    const auto num_of_total_segments = global.Width() / block_width;
    Eigen::MatrixXd A, B;
    A.resize((num_of_total_segments - num_of_ghost_segments) * block_width, (num_of_total_segments - num_of_ghost_segments) * block_width);
    A.setZero();
    B.resize((num_of_total_segments - num_of_ghost_segments) * block_width, num_of_ghost_segments * block_width);
    B.setZero();

    Vector vals;
    double tol = 1e-6;
    const auto& dense_global = global.ToDenseMatrix();
    for (auto r = 2; r < global.Height() - block_width; r++){
        dense_global->GetRow(r, vals);
        if (vals.Norml2() > tol){
            for (auto c = 0; c < vals.Size(); c++){
                if (std::abs(vals[c] - tol) > 0.0){
                    if (c < block_width){
                        B(r - block_width, c) = vals[c];
                    }
                    else if(c >= global.Width() - block_width){
                        B(r - block_width, c - A.cols()) = vals[c];
                    }
                    else{
                        A(r - block_width, c - block_width) = vals[c];
                    }
                }
            }
        }
    }
    return {A, B};
}

SGBCSolver::SGBCSolver(const SGBCProperties* sbcp, const std::pair<SGBCBoundaryInfo, SGBCBoundaryInfo>& bdrInfo) :
sbcp_(sbcp)
{

    initNodeIds(buildTargetNodeIds(sbcp->order, sbcp->num_of_segments));
    
    auto mesh = buildSGBCMesh(sbcp);
    int* partitioning = mesh.GeneratePartitioning(Mpi::WorldSize());
    auto pmesh = ParMesh(MPI_COMM_WORLD, mesh, partitioning);
    auto fec   = DG_FECollection(sbcp->order, 1, BasisType::GaussLobatto);
    auto pfes  = ParFiniteElementSpace(&pmesh, &fec);

    Model model = buildSGBCModel(mesh, partitioning, sbcp, bdrInfo);
    auto global_operator = assembleGlobalOperator(model, pfes, sbcp->order);
    std::cout << "Applying Eigen Solver on Global Operator." << std::endl;
    auto loc_ghost = buildLocalAndGhostMatrices(*global_operator, sbcp->order);
    auto es = applyEigenSolverOnLocalOperator(loc_ghost.first);
    std::cout << "Eigen Solver success." << std::endl;
    S_ = es.eigenvectors(); // S
    Sinv_ = S_.inverse(); // S-1
    F_ = Sinv_ * loc_ghost.second;
    lambda_ = es.eigenvalues(); // D = S-1 A S, flattened to eigenvalues vector.

    local_modal_values_.resize(lambda_.size());
    local_modal_values_.setZero();
    local_modal_values_old_.resize(local_modal_values_.size());
    local_modal_values_old_.setZero();
    forcing_modal_values_.resize(F_.cols());
    forcing_modal_values_.setZero();
}

SGBCNodePairInfo::SGBCNodePairInfo(const NodePair& global_pair, const size_t modal_vec_size)
{
    local_modal_values_.resize(modal_vec_size);
    local_modal_values_.setZero();
    node_pairs.g_el1 = global_pair.first;
    node_pairs.g_el2 = global_pair.second;
}

}