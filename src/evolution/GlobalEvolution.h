#pragma once

#include "solver/SourcesManager.h"
#include "solver/SolverExtension.h"
#include "components/DGOperatorFactory.h"
#include "components/Types.h"
#include <map>
#include <unordered_map>
#include <vector>

namespace maxwell {

using NodeId = int;
using NodePair = std::pair<NodeId, NodeId>;
using FieldGridFuncs = std::array<std::array<mfem::GridFunction, 3>, 2>;

class GlobalEvolution : public mfem::TimeDependentOperator {
public:
    static const int numberOfFieldComponents = 2;
    static const int numberOfMaxDimensions = 3;

    GlobalEvolution(mfem::ParFiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
    ~GlobalEvolution();

    virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;
    void ImplicitSolve(const double dt, const mfem::Vector& x, mfem::Vector& k) override;

    void commitSGBCCheckpoint(double base_time, double dt,
                              const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields);
    void finalizeSGBCStep(const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields);

    const mfem::SparseMatrix& getConstGlobalOperator() { return *globalOperator_.get(); }

    bool hasSGBC() const { return !sgbc_states_.empty(); }

private:
    void applyTFSFSourceToVector(double t_stage, int ndofs, int nbrDofs,
                                  mfem::Vector& result_vector) const;

    std::unique_ptr<mfem::SparseMatrix> globalOperator_;
    std::unique_ptr<mfem::SparseMatrix> TFSFOperator_;
    std::unique_ptr<mfem::SparseMatrix> SGBCOperator_;

    mfem::Array<int> tfsf_sub_to_parent_ids_;

    std::vector<std::unique_ptr<SGBCWrapper>> sgbcWrappers_;

    mutable std::map<GeomTag, std::vector<SGBCState>> sgbc_states_;
    std::map<GeomTag, SGBCWrapper*> sgbc_wrapper_map_;

    // Flattened list of (tag, state_index) for OpenMP-parallel sub-stepping.
    struct SGBCTask {
        GeomTag tag;
        size_t state_index;
    };
    std::vector<SGBCTask> sgbc_tasks_;
    int sgbc_omp_threads_ = 1;

    // Per-thread SGBCWrapper clones. sgbc_thread_pool_[tag][thread_id] owns
    // an independent solver so that multiple states can be advanced in parallel.
    // Thread 0 reuses the original wrapper from sgbc_wrapper_map_.
    std::map<GeomTag, std::vector<std::unique_ptr<SGBCWrapper>>> sgbc_thread_pool_;

    // Monolithic IMEX: checkpoint storage for SGBC states at start of RK4 step
    mutable std::map<GeomTag, std::vector<SGBCState>> sgbc_states_checkpoint_;
    mutable double sgbc_step_base_time_ = 0.0;
    mutable double sgbc_step_dt_ = 0.0;

    std::map<int, int> sgbc_coupling_map_;

    // Pre-computed element DOF ranges for SGBC flux gathering.
    // Maps element ID -> list of DOF indices belonging to that element.
    // Separated by pass: elem1 DOFs (pass 1) and elem2 DOFs (pass 2).
    std::vector<int> sgbc_elem1_dofs_;  // All DOFs of unique Elem1s touching SGBC faces
    std::vector<int> sgbc_elem2_dofs_;  // All DOFs of unique Elem2s touching SGBC faces

    // Cached indices of sources that are TotalField (avoids dynamic_cast per Mult)
    std::vector<int> tfsfSourceIndices_;

    mfem::ParFiniteElementSpace& fes_;
    Model& model_;
    SourcesManager& srcmngr_;
    EvolutionOptions& opts_;

    mutable std::array<mfem::ParGridFunction, 3> eOld_, hOld_;

    mutable mfem::Vector inNew_;
    mutable mfem::Vector sgbcVec_;
    mutable mfem::Vector tfsf_assembledFunc_;

    // ImplicitSolve reusable work vectors (avoid per-call allocation)
    mutable mfem::Vector implicit_inNew_;
    mutable mfem::Vector implicit_rhs_;
    mutable mfem::Vector implicit_src_;
    mutable bool implicit_work_initialized_ = false;

    // SGBC skip threshold: if both the sub-domain internal fields and the
    // global interface DOFs are below this norm, skip the sub-solve and
    // flux injection entirely for that face.
    static constexpr double sgbc_skip_threshold_ = 1e-8;

    // TFSF skip threshold: if the evaluated planewave source norm falls
    // below this value, the source has decayed and TFSF is permanently
    // skipped for the remainder of the simulation.
    static constexpr double tfsf_skip_threshold_ = 1e-8;

    // Cached dense LU factorization for small serial systems (SGBC sub-solver).
    // Activated only when n <= threshold AND nbrDofs == 0 (serial mesh).
    // This ensures it never triggers for the main 2D/3D parallel problem.
    static constexpr int dense_solve_threshold_ = 600;
    mutable mfem::DenseMatrix A_dense_;           // Dense copy of globalOperator_
    mutable mfem::DenseMatrix J_dense_;           // J = I - dt*A
    mutable std::unique_ptr<mfem::DenseMatrixInverse> J_inv_;
    mutable double cached_dt_ = -1.0;
    mutable bool dense_A_formed_ = false;
};

void load_in_to_eh_gpu(const mfem::Vector& in, 
                       std::array<mfem::ParGridFunction, 3>& e,
                       std::array<mfem::ParGridFunction, 3>& h,
                       int ndofs);
    
void load_eh_to_innew_gpu(const mfem::Vector& inOld,
                          mfem::Vector& inNew,
                          int ndofs,
                          int nbrSize);

void load_nbr_to_innew_gpu(const std::array<mfem::ParGridFunction, 3>& eOldNbr,
                   const std::array<mfem::ParGridFunction, 3>& hOldNbr,
                   mfem::Vector& inNew,
                   const int ndofs,
                   const int nbrSize);

mfem::Vector load_tfsf_into_single_vector_gpu(const FieldGridFuncs& func);

void load_tfsf_into_out_vector_gpu(const mfem::Array<int>& tfsf_sub_to_parent_ids_, 
                                   const mfem::Vector& tempTFSF,                             
                                         mfem::Vector& out,
                                   const int global_ndofs,
                                   const int tfsf_ndofs);

FieldGridFuncs eval_time_var_field_gpu(const Time time, SourcesManager& sm);

}