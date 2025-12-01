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
using NodePair = std::pair<NodeId, NodeId>;

using SGBCNodalFields = std::array<std::array<std::pair<double, double>, 3>, 2>;
using FullNodalFields = std::array<std::array<GridFunction, 3>, 2>;
using ModalValues = Eigen::VectorXcd;
using NodalValues = Eigen::VectorXd;

struct SGBCBoundaryInfo
{
    BdrCond bdrCond;
    bool isOn;
};

using SGBCBoundaries = std::pair<SGBCBoundaryInfo, SGBCBoundaryInfo>;

struct SGBCLocalNodeInfo
{
    LocalId load_el1, load_el2;
    LocalId unload_el1, unload_el2;
};

struct SGBCGlobalNodeInfo
{
    GlobalId g_el1, g_el2;
};

struct SGBCNodePairInfo
{
public:
    SGBCNodePairInfo(const NodePair& global_pair, const size_t modal_vec_size);
    
    SGBCGlobalNodeInfo node_pairs;

    void updateFullModalValues(const ModalValues& mv) { local_modal_values_ = mv; }
    void loadModalValuesInSolver(ModalValues& solver_mv);
    ModalValues getModalValues() const { return local_modal_values_; }
    
private:
    ModalValues local_modal_values_;
};

class SGBCSolver{
public:

    static std::unique_ptr<SGBCSolver> buildSGBCSolver(const SGBCProperties* sbcp); // Call to constructor call with no intBdrProperties.
    static std::unique_ptr<SGBCSolver> buildSGBCSolverWithPEC(const SGBCProperties* sbcp); // Call to constructor with PEC on both sides.

    void setFullNodalState(const FullNodalFields& in);
    FullNodalFields getFullNodalState() const;
    void setFullModalState(const ModalValues& in) { local_modal_values_ = in; }
    ModalValues getFullModalState() const { return local_modal_values_; }
    void setSGBCFieldValues(const SGBCNodalFields& in);
    SGBCNodalFields getSGBCFieldValues() const;
    void update(const Time& dt);
    size_t getLocalModalSize();

private:

    SGBCSolver(const SGBCProperties*, const SGBCBoundaries&);

    const SGBCProperties* sbcp_;
    
    SGBCLocalNodeInfo dof_pair_;
    
    Time dt_;
    
    SolverOptions opts_;

    Eigen::MatrixXcd Sinv_;
    Eigen::MatrixXcd S_;
    Eigen::MatrixXcd F_;

    ModalValues local_modal_values_, local_modal_values_old_, forcing_modal_values_;

    Eigen::VectorXcd lambda_;
    
    Fields<ParFiniteElementSpace, ParGridFunction>* global_nodal_fields_;
    
    void initNodeIds(const std::vector<NodeId>& target_ids);
    void applyLocalNodalToLocalModalTransformation(const NodalValues& in);
    void applyModalToNodalTransformation(NodalValues& out) const;
    void evol(const ModalValues& q_old, const Time& dt);

};

class SGBCWrapper
{
public:
  
  std::unique_ptr<SGBCSolver> sgbcSolver;
  std::vector<std::unique_ptr<SGBCNodePairInfo>> sgbcNodeInfos;

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);
Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnLocalOperator(const Eigen::MatrixXd& mat);

}