#include "SolverExtension.h"
#include "Solver.h"
#include "components/DGOperatorFactory.h"
#include "components/ProblemDescription.h"

namespace maxwell
{

using namespace mfem;

InteriorFaceConnectivityMaps getGlobalNodeID(const InteriorFaceConnectivityMaps& local_dof_ids, const GlobalConnectivity& global)
{
    InteriorFaceConnectivityMaps res;
    res.first.resize(local_dof_ids.first.size());
    res.second.resize(local_dof_ids.second.size());
    for (auto v{ 0 }; v < res.first.size(); v++) {
        res.first[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(local_dof_ids.first[v], local_dof_ids.second[v])));
        res.second[v] = std::distance(std::begin(global), std::find(global.begin(), global.end(), std::make_pair(local_dof_ids.second[v], local_dof_ids.first[v])));
    }
    return res;
}

void SBCSolver::findDoFPairs(Model& model, ParFiniteElementSpace& fes)
{
    // auto attMap{ mapOriginalAttributes(*fes.GetMesh()) };
    // auto fec = dynamic_cast<const DG_FECollection*>(fes.FEColl());
    // GlobalConnectivity global = assembleGlobalConnectivityMap(*fes.GetMesh(), fec);
    // auto sbc_marker = model.getMarker(BdrCond::SBC, true);
    // for (auto b = 0; b < model.getMesh().GetNBE(); b++){
    //     if (sbc_marker[model.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
    //         const FaceElementTransformations* faceTrans;
    //         fes.GetMesh()->FaceIsInterior(fes.GetMesh()->GetFaceElementTransformations(fes.GetMesh()->GetBdrElementFaceIndex(b))->ElementNo) ? faceTrans = fes.GetMesh()->GetInternalBdrFaceTransformations(b) : faceTrans = fes.GetMesh()->GetBdrFaceTransformations(b);
    //         auto twoElemSubMesh{ assembleInteriorFaceSubMesh(*fes.GetMesh(), *faceTrans, attMap) };
    //         FiniteElementSpace subFES(&twoElemSubMesh, fec);
    //         auto node_pair_global{ getGlobalNodeID(buildConnectivityForInteriorBdrFace(*faceTrans, fes, subFES), global)};
    //         for (auto p = 0; p < node_pair_global.first.size(); p++){
    //             dof_pairs_.emplace_back(node_pair_global.first[p], node_pair_global.second[p]);
    //         }
    //     }
    // }
}

GeomTagToMaterial getSBCSolverGeomTagToMaterialFromGlobal(Model& g_model)
{
    GeomTagToMaterial res;
    const auto& tag2bdr = g_model.getGeomTagToIntBoundaryCond();
    for(const auto& [tag, cond] : tag2bdr){
        if (cond == BdrCond::SBC){
            const auto sbc_material = g_model.getGeomTagToBoundaryMaterial().at(tag);
            res.emplace(1, sbc_material);
        }
    }
    return res;
}

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
    res[1] = order + 1; // start of mesh, sbc element left boundary
    res[2] = (order + 1) * (num_of_segments + 2) - (order + 1); // end of mesh, sbc element right boundary
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

Eigen::VectorXcd getEigenVectorFromOperator(const Eigen::MatrixXcd& S, int r)
{
    MFEM_ASSERT(r >= 0 && r < S.rows(), "Row is out of range.");
    return S.row(r).transpose();
}

void updateModalValues(
    const FieldComponentToFluxRows& eigvecs, 
    const Nodes& target_ids, 
    const Fields<ParFiniteElementSpace, ParGridFunction>& fields, 
    ModalValues& out)
{
    auto ndofs = fields.get(E, X).FESpace()->GetNDofs();
    const auto field_offset = 3 * ndofs;
    const auto dir_offset = ndofs;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            for (auto n = 0; n < eigvecs.at({f,d}).row_left_first.size(); n++){
                out[f * field_offset + d * dir_offset + target_ids[0]] += eigvecs.at({f,d}).row_left_first[n] * fields.get(f, d)[n];
                out[f * field_offset + d * dir_offset + target_ids[1]] += eigvecs.at({f,d}).row_left_second[n] * fields.get(f, d)[n];
                out[f * field_offset + d * dir_offset + target_ids[2]] += eigvecs.at({f,d}).row_right_first[n] * fields.get(f, d)[n];
                out[f * field_offset + d * dir_offset + target_ids[3]] += eigvecs.at({f,d}).row_right_second[n] * fields.get(f, d)[n];
            }
        }
    }
}

void updateNodalValues(const Eigen::MatrixXcd& S, const ModalValues& q, NodalValues& x)
{
    auto temp = S * q; //This should be a product giving a real vector.
    MFEM_ASSERT(temp.imag().sum() < 1e-5, "S * q product has a large imaginary module.");
    x = temp.real();
}

ModalValues evolEigenvalueSystem(const ModalValues& q, const Eigen::VectorXcd& eigvals, const Nodes& target_ids, const Time dt)
{
    ModalValues res;
    res.resize(q.size());
    res = q;
    for (auto i = 0; i < target_ids.size(); i++){
        res[target_ids[i]] = std::exp(eigvals[i] * dt) * q[i]; // + v[i], but for now lets assume no forcing sources.
    }
    return res;
}

void loadEigenVectorFromOperator(const Eigen::MatrixXcd& op, const Nodes& target_ids, const size_t ndofs, FieldComponentToFluxRows& out)
{
    const auto field_offset = 3 * ndofs;
    const auto dir_offset = ndofs;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            out[{f, d}].row_left_first   = getEigenVectorFromOperator(op, f * field_offset + d * dir_offset + target_ids[0]);
            out[{f, d}].row_left_second  = getEigenVectorFromOperator(op, f * field_offset + d * dir_offset + target_ids[1]);
            out[{f, d}].row_right_first  = getEigenVectorFromOperator(op, f * field_offset + d * dir_offset + target_ids[2]);
            out[{f, d}].row_right_second = getEigenVectorFromOperator(op, f * field_offset + d * dir_offset + target_ids[3]);
        }
    }
}

void SBCSolver::update(const Time& dt)
{
    ModalValues q_old = modal_values_;
    loadNodalValuesAtFaces();
    applyNodalToModalTransformation(); // S-1 x = q / Only flux rows belonging to the sbc interfaces
    //evol
    applyModalToNodalTransformation(); // S q = x / Full matrix-vector operation
    unloadNodalValuesAtFaces();
}

void SBCSolver::loadNodalValuesAtFaces()
{
    const auto field_offset = 3 * nodal_values_.size() / 6;
    const auto dir_offset   =     nodal_values_.size() / 6;

    for (auto f : {E, H}){ //Apply on interface nodes.
        for (auto d : {X, Y, Z}){
            nodal_values_[f * field_offset + d * dir_offset + dof_pair_.load.l_el1] = global_nodal_fields_->get(f, d)[dof_pair_.load.g_el1];
            nodal_values_[f * field_offset + d * dir_offset + dof_pair_.load.l_el2] = global_nodal_fields_->get(f, d)[dof_pair_.load.g_el2];
        }
    }

    const auto& extra_dof = dof_pair_.load.l_el1;

    for (auto f : {E, H}){ //Extend to the rest of ghost element 1.
        for (auto d : {X, Y, Z}){
            for (auto i = 0; i < extra_dof; i++){
                nodal_values_[f * field_offset + d * dir_offset + i] = nodal_values_[f * field_offset + d * dir_offset + dof_pair_.load.l_el1];
            }
        }
    }

        for (auto f : {E, H}){ //Extend to the rest of ghost element 2.
        for (auto d : {X, Y, Z}){
            for (auto i = 0; i < extra_dof; i++){
                nodal_values_[f * field_offset + d * dir_offset + dof_pair_.load.l_el2 + i] = nodal_values_[f * field_offset + d * dir_offset + dof_pair_.load.l_el2];
            }
        }
    }

}

void SBCSolver::unloadNodalValuesAtFaces()
{
    const auto field_offset = 3 * nodal_values_.size() / 6;
    const auto dir_offset   =     nodal_values_.size() / 6;

    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            global_nodal_fields_->get(f,d)[dof_pair_.unload.g_el1] = nodal_values_[f * field_offset + d * dir_offset + dof_pair_.unload.l_el1];
            global_nodal_fields_->get(f,d)[dof_pair_.unload.g_el2] = nodal_values_[f * field_offset + d * dir_offset + dof_pair_.unload.l_el2];
        }
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

void SBCSolver::initNodeIds(const std::pair<NodeId,NodeId>& el1_el2_ids)
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

SBCSolver::SBCSolver(Model& g_model, ParFiniteElementSpace& g_fes, const SBCProperties* sbcp) : 
sbcp_(sbcp)
{
    const size_t number_of_field_components = 2;
    const size_t number_of_max_dimensions = 3;
    
    target_ids_ = buildTargetNodeIds(sbcp->order, sbcp->num_of_segments);
    
    auto mesh  = Mesh::MakeCartesian1D(sbcp->num_of_segments + 2, sbcp->material_width + 2 * (sbcp->material_width / double(sbcp->num_of_segments)));
    auto pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    auto fec   = DG_FECollection(sbcp->order, 1, BasisType::GaussLobatto);
    auto pfes  = ParFiniteElementSpace(&pmesh, &fec);

    Model model(mesh, GeomTagToMaterialInfo(getSBCSolverGeomTagToMaterialFromGlobal(g_model), GeomTagToBoundaryMaterial{}), GeomTagToBoundaryInfo{});
    auto global_operator = assembleGlobalOperator(model, pfes, sbcp->order);
    auto es = applyEigenSolverOnGlobalOperator(*global_operator);
    modal_to_nodal_matrix_ = es.eigenvectors();
    const Eigen::VectorXcd D = es.eigenvalues(); // D = S-1 global_operator S, flattened to eigenvalues vector.
    
    const auto ndofs = pfes.GetNDofs();
    loadEigenVectorFromOperator(modal_to_nodal_matrix_.inverse(), target_ids_, ndofs, nodal_to_modal_rows_); // Using inverse, as q = S-1 * x

    modal_values_.resize(number_of_field_components * number_of_max_dimensions * target_ids_.size());
    modal_values_.setZero();
    nodal_values_.resize(number_of_field_components * number_of_max_dimensions * ndofs);
    nodal_values_.setZero();

    findDoFPairs(g_model, g_fes);
    
}

void SBCSolver::assignGlobalFields(Fields<ParFiniteElementSpace,ParGridFunction>* g_fields)
{
    global_nodal_fields_ = g_fields;
}

}