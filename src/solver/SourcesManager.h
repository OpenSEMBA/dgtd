#pragma once

#include "components/Sources.h"
#include "components/SubMesher.h"

#include <vector>

using FieldGridFuncs = std::array<std::array<mfem::GridFunction, 3>, 2>;

namespace maxwell {

enum class DirectSourceWaveform : int {
    Gaussian = 0,
    ModulatedGaussian = 1
};

struct DirectPlanewaveSourceData {
    int field_type = 0;
    int waveform = 0;
    double spread = 0.0;
    double mean = 0.0;
    double frequency = 0.0;
    double polarization[3] = {0.0, 0.0, 0.0};
    double propagation[3] = {0.0, 0.0, 0.0};
};

#ifdef SEMBA_DGTD_ENABLE_CUDA
void eval_tfsf_planewave_gpu(const mfem::Vector& dof_coords,
                             const mfem::Vector& tfsf_sign,
                             const mfem::Vector& source_params,
                             const mfem::Array<int>& source_meta,
                             bool apply_tfsf_sign,
                             double time,
                             int ndofs,
                             FieldGridFuncs& fields);
#endif

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
    void evalTimeVarFieldDirect(Time time, bool apply_tfsf_sign = true);
    bool hasDirectEval() const { return direct_eval_ready_; }

    /// Host-only TFSF source eval into [Ex|Ey|Ez|Hx|Hy|Hz] (6*ndofs).
    /// Does not touch cached_tfsf_fields_ / device memory — safe for MOR export
    /// without poisoning the GPU Mult path.
    void evalTimeVarFieldDirectToHost(Time time,
                                      std::vector<double>& out,
                                      bool apply_tfsf_sign = true) const;

    // Pre-computed TF/SF sign mask: +0.5 for TF DOFs, -0.5 for SF, +1 if no SF.
    const mfem::Vector& getTFSFSign() const { return tfsf_sign_; }

    // Return reference to cached field grid functions (for fast path).
    const FieldGridFuncs& getCachedTFSFFields() const { return cached_tfsf_fields_; }

    /// Find global TFSF submesh DOF whose coords match @a pos (requires initDirectPlanewaveEval).
    int findTFSFDofAtPosition(const Source::Position& pos, double tol = 1e-10) const;

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
    mfem::Vector dof_coords_;   // Flattened physical coords per DOF (3 * ndofs)
    mfem::Vector tfsf_sign_;    // TF/SF sign mask per DOF
    mfem::Vector direct_source_params_;
    mfem::Array<int> direct_source_meta_;
    int direct_source_count_ = 0;

};

}