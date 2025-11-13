#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>

namespace maxwell 
{

using namespace mfem;

using GlobalId = NodeId;
using LocalId = NodeId;

using SGBCNodalFields = std::array<std::array<std::pair<double, double>, 3>, 2>;
using FullNodalFields = std::array<std::array<GridFunction, 3>, 2>;

struct NodePairs
{
    GlobalId g_el1, g_el2;
    LocalId l_el1, l_el2;
};

struct SGBCNodeInfo
{
    NodePairs load;
    NodePairs unload;
};

using ModalValues = Eigen::VectorXcd;
using NodalValues = Eigen::VectorXd;

class SGBCSolver{
public:

    SGBCSolver(const SGBCProperties*, const std::pair<GlobalId, GlobalId>&);

    void setFullNodalState(const FullNodalFields& in);
    FullNodalFields getFullNodalState() const;
    void setSGBCFieldValues(const SGBCNodalFields& in);
    SGBCNodalFields getSGBCFieldValues() const;
    void update(const Time& dt);
    
private:
    
    const SGBCProperties* sbcp_;
    
    SGBCNodeInfo dof_pair_;
    
    Model model_;
    
    Time dt_;
    
    SolverOptions opts_;

    Eigen::MatrixXcd nodal_to_modal_matrix_;
    Eigen::MatrixXcd modal_to_nodal_matrix_;

    ModalValues modal_values_, q_old_;

    Eigen::VectorXcd eigvals_;

    std::vector<NodeId> target_ids_;
    
    Fields<ParFiniteElementSpace, ParGridFunction>* global_nodal_fields_;
    
    void initNodeIds(const std::pair<NodeId,NodeId>& el1_el2_ids);
    void applyNodalToModalTransformation(const NodalValues& in);
    void applyModalToNodalTransformation(NodalValues& out) const;
    void evol(const ModalValues& q_old, const Time& dt);

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);
Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnGlobalOperator(const SparseMatrix& mat);

}