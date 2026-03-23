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
    double rot[9];  // 3x3 rotation matrix (row-major): global→face-local frame
    
    void init(int size) {
        fields_state.SetSize(size);
        fields_state = 0.0;
        // Identity by default (correct for 1D / face-aligned cases)
        std::fill(rot, rot + 9, 0.0);
        rot[0] = rot[4] = rot[8] = 1.0;
    }
};

// Build a 3x3 rotation matrix from a face normal vector.
// Row 0 = normal (maps to sub-solver X, decoupled from 1D curl).
// Row 1 = first tangential (maps to sub-solver Y, curl-active).
// Row 2 = second tangential (maps to sub-solver Z, curl-active).
inline void buildFaceRotationMatrix(const mfem::Vector& normal, int meshDim, double rot[9])
{
    std::fill(rot, rot + 9, 0.0);
    if (meshDim <= 1) {
        rot[0] = rot[4] = rot[8] = 1.0;
        return;
    }

    double nx = normal(0);
    double ny = (normal.Size() > 1) ? normal(1) : 0.0;
    double nz = (normal.Size() > 2) ? normal(2) : 0.0;
    double len = std::sqrt(nx*nx + ny*ny + nz*nz);
    nx /= len; ny /= len; nz /= len;

    if (meshDim == 2) {
        // 2D mesh in xy-plane: z is always tangential
        rot[0] =  nx; rot[1] = ny; rot[2] = 0.0;
        rot[3] = -ny; rot[4] = nx; rot[5] = 0.0;
        rot[6] = 0.0; rot[7] = 0.0; rot[8] = 1.0;
    } else {
        // 3D: general orthonormal frame from normal
        rot[0] = nx; rot[1] = ny; rot[2] = nz;
        // First tangential: n × ref (choose ref not parallel to n)
        double rx = 0.0, ry = 0.0, rz = 0.0;
        if (std::abs(nx) < 0.9) { rx = 1.0; } else { ry = 1.0; }
        double t1x = ny*rz - nz*ry;
        double t1y = nz*rx - nx*rz;
        double t1z = nx*ry - ny*rx;
        double t1l = std::sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
        t1x /= t1l; t1y /= t1l; t1z /= t1l;
        rot[3] = t1x; rot[4] = t1y; rot[5] = t1z;
        // Second tangential: n × t1
        rot[6] = ny*t1z - nz*t1y;
        rot[7] = nz*t1x - nx*t1z;
        rot[8] = nx*t1y - ny*t1x;
    }
}

class SGBCWrapper{
public:

    static std::unique_ptr<SGBCWrapper> buildSGBCWrapper(const SGBCProperties& sbcp); 
    static std::unique_ptr<SGBCWrapper> buildSGBCWrapperWithPEC(const SGBCProperties& sbcp);

    /// Create an independent copy sharing the same SGBCProperties but with
    /// its own Solver instance. Used to enable OpenMP-parallel sub-stepping.
    std::unique_ptr<SGBCWrapper> clone() const;

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