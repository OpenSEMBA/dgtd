#pragma once

#include "components/Sources.h"
#include "components/SubMesher.h"

using FieldGridFuncs = std::array<std::array<mfem::GridFunction, 3>, 2>;

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, mfem::ParFiniteElementSpace&, Fields<ParFiniteElementSpace, ParGridFunction>& fields);

    FieldGridFuncs evalTimeVarField(const Time, FiniteElementSpace*);
    void initTFSFPreReqs(const ParMesh&, const Array<int>& marker);
    FiniteElementSpace* getTFSpace() { return tf_fes_.get(); }
    FiniteElementSpace* getSFSpace() { return sf_fes_.get(); }
    FiniteElementSpace* getGlobalTFSFSpace() { return global_tfsf_fes_.get(); }
    TotalFieldScatteredFieldSubMesher& getTFSFSubMesher() { return tfsf_submesher_; }
    void markDoFSforTForSF(FieldGridFuncs&, bool isTF);

    // Fast direct planewave evaluation (bypasses ProjectCoefficient).
    // Must be called after initTFSFPreReqs().
    void initDirectPlanewaveEval();
    void evalTimeVarFieldDirect(Time time);
    bool hasDirectEval() const { return direct_eval_ready_; }

    // Pre-computed TF/SF sign mask: +0.5 for TF DOFs, -0.5 for SF, +1 if no SF.
    const std::vector<double>& getTFSFSign() const { return tfsf_sign_; }

    // Return reference to cached field grid functions (for fast path).
    const FieldGridFuncs& getCachedTFSFFields() const { return cached_tfsf_fields_; }

    Sources sources;

private:

    void initTFSFSubMesher(const ParMesh&, const Array<int>& marker);
    void initTFSFSpaces();
    void setInitialFields(Fields<ParFiniteElementSpace, ParGridFunction>&);

    mfem::ParFiniteElementSpace& fes_;
    TotalFieldScatteredFieldSubMesher tfsf_submesher_;
    std::unique_ptr<FiniteElementSpace> tf_fes_, sf_fes_, global_tfsf_fes_;

    FieldGridFuncs cached_tfsf_fields_;
    bool cached_tfsf_fields_init_ = false;

    // Direct planewave evaluation precomputed data
    bool direct_eval_ready_ = false;
    std::vector<mfem::Vector> dof_coords_;  // Physical coords per DOF
    std::vector<double> tfsf_sign_;          // TF/SF sign mask per DOF

};

}