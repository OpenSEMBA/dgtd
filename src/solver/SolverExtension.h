#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>
#include <memory>

namespace maxwell 
{

class Solver;

using GlobalId = NodeId;
using LocalId = NodeId;
using GhostId = NodeId;
using NodePair = std::pair<NodeId, NodeId>;
using ElementPair = std::pair<ElementId, ElementId>;

using SGBCForcingFields = std::array<std::array<std::pair<double, double>, 3>, 2>;
using SGBCNodalFields = SGBCForcingFields;
using FullNodalFields = std::array<std::array<mfem::GridFunction, 3>, 2>;
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

struct SGBCState {
    NodePair global_pair;       
    ElementPair element_pair;   
    mfem::Vector fields_state;  
    
    void init(int size) {
        fields_state.SetSize(size);
        fields_state = 0.0;
    }
};

class SGBCWrapper{
public:

    static std::unique_ptr<SGBCWrapper> buildSGBCWrapper(const SGBCProperties& sbcp); 
    static std::unique_ptr<SGBCWrapper> buildSGBCWrapperWithPEC(const SGBCProperties& sbcp);

    // [MODIFIED] Now takes a specific state context
    void updateFieldsWithGlobal(const std::array<mfem::ParGridFunction, 3>& e, 
                                const std::array<mfem::ParGridFunction, 3>& h, 
                                const SGBCState& context);

    void setAllSolverFields(const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields);
    
    // [MODIFIED] Now reads from state directly
    void getSGBCFields(const Array<int>& sub_to_global, const SGBCState& context, FieldGridFuncs& out);
    
    const SGBCProperties& getProperties() const { return sbcp_; }

    void solve(const Time t, const Time dt);
    void setOldTime(const Time t) { old_t_ = t; }
    const Time getOldTime() const { return old_t_; }

    // [ADDED] Context Switching
    void loadState(const SGBCState& state);
    void saveState(SGBCState& state);
    int getStateSize() const;

    ~SGBCWrapper();

private:

    SGBCWrapper(const SGBCProperties&, const SGBCBoundaries&);

    const SGBCProperties& sbcp_;
    
    std::unique_ptr<Solver> solver_;

    SGBCLocalNodeInfo dof_pair_;
    FluxNodeInfo flux_nodes_;

    Time old_t_ = 0.0;
    
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>* global_fields_;
    std::unique_ptr<Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>> sgbc_solver_fields_;
    
    void initNodeIds(const std::vector<NodeId>& target_ids);

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);

}