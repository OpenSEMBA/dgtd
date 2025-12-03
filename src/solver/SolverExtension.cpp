#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

namespace maxwell
{

using namespace mfem;

const auto num_of_field_components = 2;
const auto num_of_max_dimensions = 3;
const auto num_of_field_blocks = num_of_field_components * num_of_max_dimensions;
const auto num_of_ghost_segments_per_field_comp = 2;

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

std::vector<NodeId> buildTargetNodeIds(const size_t order, const size_t num_of_segments)
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

    int dofs_per_element = sbcp_->order + 1;
    int ndofs = sbcp_->num_of_segments * dofs_per_element;
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

void SGBCSolver::setForcingNodalToModalFieldValues(const SGBCForcingFields& in)
{
    NodalValues forcing_nodal_values(F_.cols());
    const auto num_of_entries_per_block = forcing_nodal_values.size() / num_of_field_blocks;
    const auto ghost_block_size = num_of_entries_per_block / num_of_ghost_segments_per_field_comp;
    
    const auto field_offset = num_of_entries_per_block * num_of_max_dimensions;
    const auto dir_offset = num_of_entries_per_block;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            for (auto v = 0; v < ghost_block_size; v++){
                forcing_nodal_values[f * field_offset + d * dir_offset                    + v] = in.at(f).at(d).first; // First Ghost Element - Same value for all nodes.
                forcing_nodal_values[f * field_offset + d * dir_offset + ghost_block_size + v] = in.at(f).at(d).second; // Second Ghost Element - Same value for all nodes.
            }
        }
    }

    forcing_modal_values_ = F_ * forcing_nodal_values;

}

void SGBCSolver::update(const Time& dt)
{
    local_modal_values_old_ = local_modal_values_;
    evol(local_modal_values_old_, dt);
}

void SGBCSolver::evol(const ModalValues& local_modal_values_old_, const Time& dt)
{
    // Eigen::VectorXcd forcing = forcing_modal_values_;
    for (auto i = 0; i < local_modal_values_.size(); i++){
        local_modal_values_[i] = std::exp(lambda_[i] * dt) * local_modal_values_old_[i] + forcing_modal_values_[i];
    }
}

void SGBCSolver::applyLocalNodalToLocalModalTransformation(const NodalValues& in)
{
    local_modal_values_ = Sinv_ * in;
}

void SGBCSolver::applyModalToNodalTransformation(NodalValues& out) const
{
    const Eigen::VectorXcd temp = S_ * local_modal_values_;
    for (auto v = 0; v < temp.size(); v++){
        if (temp[v].imag() > 1e-5){
            throw std::runtime_error("Imaginary part over tolerance in ModalToNodalTransformation.");
        }
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
    Material vacuum = buildVacuumMaterial();
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

bool matches_set(const std::unordered_set<int>& baseSet, int n, int offset) {
    int r = n % offset;
    return baseSet.find(r) != baseSet.end();
}

std::vector<int> getGlobalToGhostElementColumns(const int global_size, const int dofs_per_element, bool isLeft)
{
    const auto num_of_field_components = 2;
    const auto num_of_max_dimensions = 3;
    const auto num_of_ghost_segments = 2;
    const auto field_blocks = num_of_field_components * num_of_max_dimensions;
    const auto field_block_size = global_size / field_blocks;

    std::vector<int> ghosts;
    for (auto d = 0; d < dofs_per_element; d++){
        if (isLeft){
            ghosts.push_back(d);
        }
        else{
            ghosts.push_back(field_block_size - d - 1);
        }
    }
    std::unordered_set<int> ghostSet(ghosts.begin(), ghosts.end());
    
    std::vector<int> res;
    for(auto c = 0; c < global_size; c++){
        if(matches_set(ghostSet, c, field_block_size)){
            res.push_back(c);
        }
    }
    return res;
}

std::vector<int> getGlobalToRealElementColumns(const std::vector<int>& left_ghost, const std::vector<int>& right_ghost, const int global_size)
{
    std::vector<int> res;
    std::unordered_set<int> ghost_set;

    ghost_set.insert(left_ghost.begin(), left_ghost.end());
    ghost_set.insert(right_ghost.begin(), right_ghost.end());

    for (int i = 0; i < global_size; ++i) {
        if (ghost_set.find(i) == ghost_set.end()) {
            res.push_back(i);
        }
    }
    return res;
}

std::unordered_map<int, size_t> buildColumnsForOperatorB(const std::vector<int>& left_ghost_cols, const std::vector<int>& right_ghost_cols){
    std::unordered_map<int, size_t> res;
    res.reserve(left_ghost_cols.size() + right_ghost_cols.size());

    size_t total = left_ghost_cols.size() + right_ghost_cols.size();
    size_t i = 0, j = 0;
    size_t next_b_col = 0;

    for (size_t k = 0; k < total; ++k) {
        bool take_left =
            (j >= right_ghost_cols.size()) ||
            (i < left_ghost_cols.size() && left_ghost_cols[i] < right_ghost_cols[j]);

        if (take_left) {
            res[left_ghost_cols[i]] = next_b_col++;
            ++i;
        } else {
            res[right_ghost_cols[j]] = next_b_col++;
            ++j;
        }
    }
    return res;
}

std::unordered_map<int, size_t> buildColumnsForOperatorA(const std::vector<int>& real_cols){
    std::unordered_map<int, size_t> res;
    res.reserve(real_cols.size());
    for (size_t k = 0; k < real_cols.size(); ++k) {
        res[real_cols[k]] = k;
    }
    return res;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> buildLocalAndGhostMatrices(const SparseMatrix& global, const int order)
{
    const auto field_blocks = num_of_field_components * num_of_max_dimensions;
    const auto field_block_size = global.Width() / field_blocks;
    const auto num_of_total_ghost_segments = 2 * field_blocks;
    const auto dofs_per_element = order + 1;
    const auto num_of_total_segments = global.Width() / dofs_per_element;

    auto left_ghost_cols = getGlobalToGhostElementColumns(global.Width(), dofs_per_element, true);
    auto right_ghost_cols = getGlobalToGhostElementColumns(global.Width(), dofs_per_element, false);
    auto real_cols = getGlobalToRealElementColumns(left_ghost_cols, right_ghost_cols, global.Width());

    auto b_col_map = buildColumnsForOperatorB(left_ghost_cols, right_ghost_cols);
    auto a_col_map = buildColumnsForOperatorA(real_cols);
    
    int row_offset = 0;
    Vector vals;
    Eigen::MatrixXd A, B;
    A.resize((num_of_total_segments - num_of_total_ghost_segments) * dofs_per_element, (num_of_total_segments - num_of_total_ghost_segments) * dofs_per_element);
    A.setZero();
    B.resize((num_of_total_segments - num_of_total_ghost_segments) * dofs_per_element, num_of_total_ghost_segments * dofs_per_element);
    B.setZero();
    const auto& dense_global = global.ToDenseMatrix();
    for (auto r = 0; r < global.Height(); ++r) {
        if (std::find(real_cols.begin(), real_cols.end(), r) == real_cols.end()) {
            ++row_offset;
            continue;
        }
        dense_global->GetRow(r, vals);
        for (auto c = 0; c < global.Width(); ++c) {
            auto itB = b_col_map.find(c);
            if (itB != b_col_map.end()) {
                B(r - row_offset, itB->second) = vals[c];
                continue;
            }
            auto itA = a_col_map.find(c);
            if (itA != a_col_map.end()) {
                A(r - row_offset, itA->second) = vals[c];
            }
        }
    }

    return {A, B};
}

SGBCNodalFields SGBCSolver::getSGBCNodalFields() const
{
    const auto dofs_per_element = sbcp_->order + 1;
    const auto local_nodal_values = getFullNodalState();
    SGBCNodalFields res;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            res.at(f).at(d).first = local_nodal_values.at(f).at(d)[dofs_per_element];
            res.at(f).at(d).second = local_nodal_values.at(f).at(d)[local_nodal_values.at(f).at(d).Size() - dofs_per_element - 1];
        }
    }
    return res;
}

SGBCSolver::SGBCSolver(const SGBCProperties* sbcp, const std::pair<SGBCBoundaryInfo, SGBCBoundaryInfo>& bdrInfo) :
sbcp_(sbcp)
{ 
    
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

    const auto forcing_nodes_field_block_size = F_.cols() / num_of_field_blocks;
    const auto local_nodes_field_block_size = local_modal_values_.size() / num_of_field_blocks;
    flux_nodes_.ghost_element_left_node = sbcp_->order;
    flux_nodes_.ghost_element_right_node = forcing_nodes_field_block_size - 1 - sbcp_->order;
    flux_nodes_.local_element_left_node = 0;
    flux_nodes_.local_element_right_node = local_nodes_field_block_size - 1 - sbcp_->order;

    // initNodeIds(buildTargetNodeIds(sbcp->order, sbcp->num_of_segments));
}

SGBCNodePairInfo::SGBCNodePairInfo(const NodePair& global_pair, const size_t modal_vec_size)
{
    local_modal_values_.resize(modal_vec_size);
    local_modal_values_.setZero();
    node_pairs.g_el1 = global_pair.first;
    node_pairs.g_el2 = global_pair.second;
}

}