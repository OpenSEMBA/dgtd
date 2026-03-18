#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>
#include <memory>

namespace maxwell 
{

class Solver;

using SGBCHelperFields = std::array<std::array<mfem::Vector, 3>, 2>;
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
    void updateFieldsWithGlobal(const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields, 
                                const SGBCState& context);

    // Overload that reads ghost-zone values from a raw stage vector (for monolithic IMEX)
    void updateFieldsWithGlobalVector(const mfem::Vector& in, int ndofs, const SGBCState& context);

    void setAllSolverFields(const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields);
    
    // [MODIFIED] Now reads from state directly
    void getSGBCFields(const Array<int>& sub_to_global, const SGBCState& context, SGBCHelperFields& out);

    // Write SGBC sub-solver boundary values directly into a global-layout vector
    void fillGlobalSGBCVec(const SGBCState& context, mfem::Vector& vec, int blockSize);
    
    const SGBCProperties& getProperties() const { return sbcp_; }

    void solve(const Time t, const Time dt);
    void setOldTime(const Time t) { old_t_ = t; }
    const Time getOldTime() const { return old_t_; }
    double getRecommendedDt() const { return recommended_dt_; }

    // [ADDED] Context Switching
    void loadState(const SGBCState& state);
    void saveState(SGBCState& state);
    int getStateSize() const;

    // Interface DOF indices for the sub-solver's internal field layout
    int getLocalFieldSize() const;
    int getLeftInterfaceIndex() const;
    int getRightInterfaceIndex() const;

    ~SGBCWrapper();

private:

    SGBCWrapper(const SGBCProperties&, const SGBCBoundaries&);

    const SGBCProperties& sbcp_;
    
    std::unique_ptr<Solver> solver_;

    SGBCLocalNodeInfo dof_pair_;
    FluxNodeInfo flux_nodes_;

    Time old_t_ = 0.0;
    double recommended_dt_ = 0.0;
    int n_ghost_elements_ = 1;

    bool temporal_warning_printed_ = false;
    
    const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>* global_fields_;
    std::unique_ptr<Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>> sgbc_solver_fields_;

};

}