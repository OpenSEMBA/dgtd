#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>

namespace maxwell 
{

using namespace mfem;

GeomTagToMaterial getSBCSolverGeomTagToMaterialFromGlobal(Model& global_model);

using GlobalId = NodeId;
using LocalId = NodeId;

struct FluxRows
{
    Eigen::VectorXcd row_left_first;
    Eigen::VectorXcd row_left_second;
    Eigen::VectorXcd row_right_first;
    Eigen::VectorXcd row_right_second;
};

struct NodePairs
{
    GlobalId g_el1, g_el2;
    LocalId l_el1, l_el2;
};

struct SBCNodeInfo
{
    NodePairs load;
    NodePairs unload;
};

using FieldComponentToFluxRows = std::map<std::pair<FieldType,Direction>,FluxRows>;
using ModalValues = Eigen::VectorXcd;
using NodalValues = Eigen::VectorXd;

class SBCSolver{
public:

    SBCSolver(const SBCProperties*, const std::pair<GlobalId, GlobalId>&);

    void update(const Time& dt);
    
private:
    
    const SBCProperties* sbcp_;
    
    SBCNodeInfo dof_pair_;
    
    Model model_;
    
    Time dt_;
    
    SolverOptions opts_;

    FieldComponentToFluxRows nodal_to_modal_rows_;
    Eigen::MatrixXcd modal_to_nodal_matrix_;

    ModalValues modal_values_, q_old_;
    NodalValues nodal_values_;

    Eigen::VectorXcd eigvals_;

    std::vector<NodeId> target_ids_;
    
    Fields<ParFiniteElementSpace, ParGridFunction>* global_nodal_fields_;
    
    void initNodeIds(const std::pair<NodeId,NodeId>& el1_el2_ids);
    void applyNodalToModalTransformation();
    void applyModalToNodalTransformation();
    void evol(const ModalValues& q_old, const Time& dt);

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);
Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnGlobalOperator(const SparseMatrix& mat);

}