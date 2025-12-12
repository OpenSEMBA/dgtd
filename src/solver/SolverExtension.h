#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>

namespace maxwell 
{

using GlobalId = NodeId;
using LocalId = NodeId;
using GhostId = NodeId;
using NodePair = std::pair<NodeId, NodeId>;

using SGBCForcingFields = std::array<std::array<std::pair<double, double>, 3>, 2>;
using SGBCNodalFields = SGBCForcingFields;
using FullNodalFields = std::array<std::array<GridFunction, 3>, 2>;
using ModalValues = Eigen::VectorXcd;
using NodalValues = Eigen::VectorXd;

struct SGBCBoundaryInfo
{
    BdrCond bdrCond = BdrCond::SMA;
    bool isOn = false;
};

using SGBCBoundaries = std::pair<SGBCBoundaryInfo, SGBCBoundaryInfo>;

struct FluxNodeInfo
{
    LocalId local_element_left_node, local_element_right_node;  
    GhostId ghost_element_left_node, ghost_element_right_node;
};

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
    SGBCNodePairInfo(const NodePair& global_pair);
    
    SGBCGlobalNodeInfo node_pairs;

};

class SGBCWrapper{
public:

    static std::unique_ptr<SGBCWrapper> buildSGBCWrapper(const SGBCProperties& sbcp); // Call to constructor call with no intBdrProperties.
    static std::unique_ptr<SGBCWrapper> buildSGBCWrapperWithPEC(const SGBCProperties& sbcp); // Call to constructor with PEC on both real/ghost interfaces.

    ~SGBCWrapper();

private:

    SGBCWrapper(const SGBCProperties&, const SGBCBoundaries&);

    const SGBCProperties& sbcp_;

    std::unique_ptr<Solver> solver_;

    SGBCLocalNodeInfo dof_pair_;
    FluxNodeInfo flux_nodes_;
    
    Time dt_;
    
    Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>* global_nodal_fields_;
    
    void initNodeIds(const std::vector<NodeId>& target_ids);

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);

}