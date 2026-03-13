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
    
    virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;
    void ImplicitSolve(const double dt, const mfem::Vector& x, mfem::Vector& k) override;

    void advanceSGBCs(double time, double dt,
                      const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields);

    const mfem::SparseMatrix& getConstGlobalOperator() { return *globalOperator_.get(); }

    bool hasSGBC() const { return !sgbc_states_.empty(); }

private:
    void applyTFSFSourceToVector(double t_stage, int ndofs, int ndofs_tfsf,
                                  mfem::Vector& result_vector, bool check_zero = false) const;

    std::unique_ptr<mfem::SparseMatrix> globalOperator_;
    std::unique_ptr<mfem::SparseMatrix> TFSFOperator_;
    std::unique_ptr<mfem::SparseMatrix> SGBCOperator_;
    int SGBCndofs_;

    mfem::Array<int> tfsf_sub_to_parent_ids_;
    mfem::Array<int> sgbc_sub_to_parent_ids_;

    std::vector<std::unique_ptr<SGBCWrapper>> sgbcWrappers_;

    std::map<GeomTag, std::vector<SGBCState>> sgbc_states_;
    std::map<GeomTag, SGBCWrapper*> sgbc_wrapper_map_;

    std::map<int, int> sgbc_coupling_map_;

    // O(1) lookup: global DOF id -> submesh index (replaces linear Array::Find)
    std::unordered_map<int, int> sgbc_parent_to_sub_;

    // Cached indices of sources that are TotalField (avoids dynamic_cast per Mult)
    std::vector<int> tfsfSourceIndices_;

    // Fast-exit flag: set to true once TFSF cutoff time is reached
    mutable bool tfsf_cutoff_reached_ = false;

    // Precomputed SGBC scatter arrays: sgbc_scatter_out_[i] = parent DOF index for submesh DOF i
    std::vector<int> sgbc_scatter_parent_idx_;

    mfem::ParFiniteElementSpace& fes_;
    Model& model_;
    SourcesManager& srcmngr_;
    EvolutionOptions& opts_;

    mutable std::array<mfem::ParGridFunction, 3> eOld_, hOld_;
    mutable SGBCHelperFields sgbc_helper_fields_;
    mutable int last_sgbc_helper_size_ = -1;

    mutable mfem::Vector inNew_;
    mutable mfem::Vector sgbcVec_;
    mutable mfem::Vector tempSGBC_;
    mutable mfem::Vector tfsf_assembledFunc_;
    mutable mfem::Vector tfsf_tempVec_;

    // ImplicitSolve reusable work vectors (avoid per-call allocation)
    mutable mfem::Vector implicit_inNew_;
    mutable mfem::Vector implicit_rhs_;
    mutable mfem::Vector implicit_src_;
    mutable bool implicit_work_initialized_ = false;
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