#include "SourcesManager.h"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <vector>

namespace maxwell {

using namespace mfem;

namespace {

constexpr int kDirectSourceMetaStride = 2;
constexpr int kDirectSourceParamStride = 9;
constexpr double kPi = 3.14159265358979323846;

int metaIndex(int source_idx, int offset)
{
    return source_idx * kDirectSourceMetaStride + offset;
}

int paramIndex(int source_idx, int offset)
{
    return source_idx * kDirectSourceParamStride + offset;
}

bool isClosedBasis(const FiniteElementSpace& fes)
{
    const auto* dg_fec = dynamic_cast<const mfem::DG_FECollection*>(fes.FEColl());
    return dg_fec != nullptr &&
           dg_fec->GetBasisType() == mfem::BasisType::GaussLobatto;
}

bool extractDirectPlanewaveSource(const TotalField& tf, DirectPlanewaveSourceData& out)
{
    const auto* planewave = dynamic_cast<const Planewave*>(tf.function());
    if (planewave == nullptr) {
        return false;
    }

    const auto* gaussian = dynamic_cast<const Gaussian*>(planewave->function());
    const auto* modulated = dynamic_cast<const ModulatedGaussian*>(planewave->function());
    if (gaussian == nullptr && modulated == nullptr) {
        return false;
    }

    out.field_type = static_cast<int>(planewave->fieldType());
    out.waveform = static_cast<int>(gaussian != nullptr
        ? DirectSourceWaveform::Gaussian
        : DirectSourceWaveform::ModulatedGaussian);
    out.spread = gaussian != nullptr ? gaussian->spread() : modulated->spread();
    out.mean = gaussian != nullptr ? gaussian->mean()[X] : modulated->mean()[X];
    out.frequency = modulated != nullptr ? modulated->frequency() : 0.0;
    for (int d = 0; d < 3; ++d) {
        out.polarization[d] = planewave->polarization()[d];
        out.propagation[d] = planewave->propagation()[d];
    }
    return true;
}

double evalDirectPlanewaveComponent(const DirectPlanewaveSourceData& src,
                                    double x, double y, double z,
                                    double time, int field_type, int direction)
{
    const double phase_delay =
        (x * src.propagation[0] + y * src.propagation[1] + z * src.propagation[2]) /
        physicalConstants::speedOfLight;
    const double arg = phase_delay - time - src.mean;
    const double envelope =
        std::exp(-(arg * arg) / (2.0 * src.spread * src.spread));
    const double waveform =
        (src.waveform == static_cast<int>(DirectSourceWaveform::ModulatedGaussian))
            ? envelope * std::cos(2.0 * kPi * src.frequency * arg)
            : envelope;

    double cross[3];
    if (src.field_type == static_cast<int>(E)) {
        cross[0] = src.propagation[1] * src.polarization[2] - src.propagation[2] * src.polarization[1];
        cross[1] = src.propagation[2] * src.polarization[0] - src.propagation[0] * src.polarization[2];
        cross[2] = src.propagation[0] * src.polarization[1] - src.propagation[1] * src.polarization[0];
        return waveform * (field_type == static_cast<int>(E)
            ? src.polarization[direction]
            : cross[direction]);
    }

    cross[0] = src.polarization[1] * src.propagation[2] - src.polarization[2] * src.propagation[1];
    cross[1] = src.polarization[2] * src.propagation[0] - src.polarization[0] * src.propagation[2];
    cross[2] = src.polarization[0] * src.propagation[1] - src.polarization[1] * src.propagation[0];
    return waveform * (field_type == static_cast<int>(H)
        ? src.polarization[direction]
        : cross[direction]);
}

} // namespace

SourcesManager::SourcesManager(const Sources& srcs, mfem::ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>& fields) :
    sources{ srcs }, 
    fes_{ fes }
{
    setInitialFields(fields);
}

void SourcesManager::setInitialFields(Fields<ParFiniteElementSpace, ParGridFunction>& fields)
{
    for (const auto& source : sources) {
        auto src{ dynamic_cast<InitialField*>(source.get()) };
        if (src == nullptr) {
            continue;
        }
        for (FieldType ft: {E, H}) {
            for (auto x : { X, Y, Z }) {
                auto f = [src, ft, x](const Source::Position& pos, Source::Time t) {
                    return src->eval(pos, t, ft, x);
                };
                FunctionCoefficient fc(f);
                GridFunction gf(fields.get(ft, x).FESpace());
                gf.UseDevice(true);
                gf.ProjectCoefficient(fc);
                fields.get(ft, x) += gf;
            }
        }
    }
}

FieldGridFuncs SourcesManager::evalTimeVarField(const Time time, FiniteElementSpace* fes)
{
    if (!cached_tfsf_fields_init_) {
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                cached_tfsf_fields_[ft][d].UseDevice(true);
                cached_tfsf_fields_[ft][d].SetSpace(fes);
            }
        }
        cached_tfsf_fields_init_ = true;
    }

    for (auto ft : { E, H }) {
        for (auto d : { X, Y, Z }) {
            cached_tfsf_fields_[ft][d] = 0.0;
        }
    }

    for (const auto& source : sources) {
        auto tf = dynamic_cast<TotalField*>(source.get());
        if (tf == nullptr) {
            continue;
        }
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                auto f = [tf, ft, d](const Source::Position& pos, Source::Time t) {
                    return tf->eval(pos, t, ft, d);
                };
                FunctionCoefficient func(f);
                func.SetTime(time);
                cached_tfsf_fields_[ft][d].ProjectCoefficient(func);
            }
        }
    }
    return cached_tfsf_fields_;
}

void SourcesManager::markDoFSforTForSF(FieldGridFuncs& gfs, bool isTF)
{
    auto global_tfsf_map = tfsf_submesher_.getGlobalTFSFSubMesh()->GetParentElementIDMap();
    Array<int> secondary_map;
    switch (isTF) {
    case true:
        if (tfsf_submesher_.getSFSubMesh() != NULL) {
            secondary_map = tfsf_submesher_.getSFSubMesh()->GetParentElementIDMap();
        }
        break;
    case false:
        if (tfsf_submesher_.getTFSubMesh() != NULL) {
            secondary_map = tfsf_submesher_.getTFSubMesh()->GetParentElementIDMap();
        }
        break;
    }

    for (int e = 0; e < secondary_map.Size(); e++) {
        Array<int> dofs;
        global_tfsf_fes_->GetElementDofs(global_tfsf_map.Find(secondary_map[e]), dofs);
        for (int i = 0; i < dofs.Size(); i++) {
            for (auto f : { E, H }) {
                for (auto d{ X }; d <= Z; d++) {
                    gfs[f][d].UseDevice(true);
                    gfs[f][d][dofs[i]] = 0.0;
                }
            }
        }
    }
}

void SourcesManager::initTFSFPreReqs(const ParMesh& m, const Array<int>& marker)
{
    initTFSFSubMesher(m, marker);
    initTFSFSpaces();
}

void SourcesManager::initTFSFSubMesher(const ParMesh& m, const Array<int>& marker)
{
    auto sm = TotalFieldScatteredFieldSubMesher(m, marker);
    tfsf_submesher_ = std::move(sm);
}

void SourcesManager::initTFSFSpaces()
{
    if (tfsf_submesher_.getTFSubMesh() != NULL) {
        tf_fes_ = std::make_unique<FiniteElementSpace>(tfsf_submesher_.getTFSubMesh(), fes_.FEColl());
    }
    if (tfsf_submesher_.getSFSubMesh() != NULL) {
        sf_fes_ = std::make_unique<FiniteElementSpace>(tfsf_submesher_.getSFSubMesh(), fes_.FEColl());
    }
    global_tfsf_fes_ = std::make_unique<FiniteElementSpace>(tfsf_submesher_.getGlobalTFSFSubMesh(), fes_.FEColl());
}

void SourcesManager::initDirectPlanewaveEval()
{
    direct_eval_ready_ = false;
    direct_source_count_ = 0;
    direct_source_params_.SetSize(0);
    direct_source_meta_.SetSize(0);

    if (!global_tfsf_fes_) return;
    auto* fes = global_tfsf_fes_.get();
    const int ndofs = fes->GetNDofs();
    if (!isClosedBasis(*fes)) return;

    std::vector<DirectPlanewaveSourceData> direct_sources;
    direct_sources.reserve(sources.size());
    for (const auto& source : sources) {
        auto* tf = dynamic_cast<TotalField*>(source.get());
        if (tf == nullptr) {
            continue;
        }
        DirectPlanewaveSourceData src_data;
        if (!extractDirectPlanewaveSource(*tf, src_data)) {
            return;
        }
        direct_sources.push_back(src_data);
    }
    if (direct_sources.empty()) return;

    // 1) Pre-compute physical coordinates for each DOF
    dof_coords_.SetSize(3 * ndofs);
    dof_coords_ = 0.0;
    for (int el = 0; el < fes->GetNE(); ++el) {
        const auto* fe = fes->GetFE(el);
        auto* T = fes->GetElementTransformation(el);
        const auto& ir = fe->GetNodes();
        mfem::Array<int> dofs;
        fes->GetElementDofs(el, dofs);
        for (int i = 0; i < dofs.Size(); ++i) {
            const auto& ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);
            mfem::Vector coords;
            T->Transform(ip, coords);
            const int dof = dofs[i];
            for (int c = 0; c < coords.Size(); ++c) {
                dof_coords_[3 * dof + c] = coords[c];
            }
        }
    }

    // 2) Pre-compute TF/SF sign mask
    tfsf_sign_.SetSize(ndofs);
    tfsf_sign_ = 1.0;  // Default: no SF submesh, scale = 1.0
    if (tfsf_submesher_.getSFSubMesh() != NULL) {
        // Build set of SF DOF indices (on the global TFSF submesh)
        auto global_tfsf_map = tfsf_submesher_.getGlobalTFSFSubMesh()->GetParentElementIDMap();

        // Build reverse lookup: parent_elem_id -> global_tfsf_element_index
        std::unordered_map<int, int> parent_to_global;
        for (int i = 0; i < global_tfsf_map.Size(); ++i) {
            parent_to_global[global_tfsf_map[i]] = i;
        }

        // All DOFs start as TF (+0.5)
        tfsf_sign_ = 0.5;

        // Mark SF DOFs as -0.5
        auto sf_map = tfsf_submesher_.getSFSubMesh()->GetParentElementIDMap();
        std::unordered_set<int> sf_dof_set;
        for (int e = 0; e < sf_map.Size(); ++e) {
            auto it = parent_to_global.find(sf_map[e]);
            if (it == parent_to_global.end()) continue;
            mfem::Array<int> dofs;
            fes->GetElementDofs(it->second, dofs);
            for (int i = 0; i < dofs.Size(); ++i) {
                sf_dof_set.insert(dofs[i]);
            }
        }
        for (int dof : sf_dof_set) {
            tfsf_sign_[dof] = -0.5;
        }
    }

    direct_source_count_ = static_cast<int>(direct_sources.size());
    direct_source_meta_.SetSize(direct_source_count_ * kDirectSourceMetaStride);
    direct_source_params_.SetSize(direct_source_count_ * kDirectSourceParamStride);
    for (int s = 0; s < direct_source_count_; ++s) {
        const auto& src = direct_sources[s];
        direct_source_meta_[metaIndex(s, 0)] = src.field_type;
        direct_source_meta_[metaIndex(s, 1)] = src.waveform;
        direct_source_params_[paramIndex(s, 0)] = src.spread;
        direct_source_params_[paramIndex(s, 1)] = src.mean;
        direct_source_params_[paramIndex(s, 2)] = src.frequency;
        for (int d = 0; d < 3; ++d) {
            direct_source_params_[paramIndex(s, 3 + d)] = src.polarization[d];
            direct_source_params_[paramIndex(s, 6 + d)] = src.propagation[d];
        }
    }

    // 3) Initialize cached grid functions if not already done
    if (!cached_tfsf_fields_init_) {
        for (auto ft : { E, H }) {
            for (auto d : { X, Y, Z }) {
                cached_tfsf_fields_[ft][d].UseDevice(true);
                cached_tfsf_fields_[ft][d].SetSpace(fes);
            }
        }
        cached_tfsf_fields_init_ = true;
    }

    direct_eval_ready_ = true;
}

void SourcesManager::evalTimeVarFieldDirect(Time time, bool apply_tfsf_sign)
{
    auto* fes = global_tfsf_fes_.get();
    const int ndofs = fes->GetNDofs();
    const bool on_cuda =
#ifdef SEMBA_DGTD_ENABLE_CUDA
        mfem::Device::Allows(mfem::Backend::CUDA);
#else
        false;
#endif
    if (on_cuda) {
#ifdef SEMBA_DGTD_ENABLE_CUDA
        eval_tfsf_planewave_gpu(dof_coords_,
                                tfsf_sign_,
                                direct_source_params_,
                                direct_source_meta_,
                                apply_tfsf_sign,
                                time,
                                ndofs,
                                cached_tfsf_fields_);
#endif
        return;
    }

    for (auto ft : { E, H }) {
        for (auto d : { X, Y, Z }) {
            cached_tfsf_fields_[ft][d] = 0.0;
        }
    }

    const double* coords = dof_coords_.HostRead();
    const double* sign = tfsf_sign_.HostRead();
    const double* params = direct_source_params_.HostRead();
    const int* meta = direct_source_meta_.HostRead();
    for (int i = 0; i < ndofs; ++i) {
        const double x = coords[3 * i + 0];
        const double y = coords[3 * i + 1];
        const double z = coords[3 * i + 2];
        const double scale = apply_tfsf_sign ? sign[i] : 1.0;
        for (int s = 0; s < direct_source_count_; ++s) {
            DirectPlanewaveSourceData src;
            src.field_type = meta[metaIndex(s, 0)];
            src.waveform = meta[metaIndex(s, 1)];
            src.spread = params[paramIndex(s, 0)];
            src.mean = params[paramIndex(s, 1)];
            src.frequency = params[paramIndex(s, 2)];
            for (int d = 0; d < 3; ++d) {
                src.polarization[d] = params[paramIndex(s, 3 + d)];
                src.propagation[d] = params[paramIndex(s, 6 + d)];
            }
            for (int ft : { E, H }) {
                for (int d : { X, Y, Z }) {
                    cached_tfsf_fields_[ft][d][i] +=
                        scale * evalDirectPlanewaveComponent(src, x, y, z, time, ft, d);
                }
            }
        }
    }
}

void SourcesManager::evalTimeVarFieldDirectToHost(Time time,
                                                  std::vector<double>& out,
                                                  bool apply_tfsf_sign) const
{
    MFEM_VERIFY(direct_eval_ready_, "initDirectPlanewaveEval() must be called first");
    auto* fes = global_tfsf_fes_.get();
    const int ndofs = fes->GetNDofs();
    out.assign(static_cast<std::size_t>(6 * ndofs), 0.0);

    const double* coords = dof_coords_.HostRead();
    const double* sign = tfsf_sign_.HostRead();
    const double* params = direct_source_params_.HostRead();
    const int* meta = direct_source_meta_.HostRead();

    for (int i = 0; i < ndofs; ++i) {
        const double x = coords[3 * i + 0];
        const double y = coords[3 * i + 1];
        const double z = coords[3 * i + 2];
        const double scale = apply_tfsf_sign ? sign[i] : 1.0;
        for (int s = 0; s < direct_source_count_; ++s) {
            DirectPlanewaveSourceData src;
            src.field_type = meta[metaIndex(s, 0)];
            src.waveform = meta[metaIndex(s, 1)];
            src.spread = params[paramIndex(s, 0)];
            src.mean = params[paramIndex(s, 1)];
            src.frequency = params[paramIndex(s, 2)];
            for (int d = 0; d < 3; ++d) {
                src.polarization[d] = params[paramIndex(s, 3 + d)];
                src.propagation[d] = params[paramIndex(s, 6 + d)];
            }
            for (int ft : { E, H }) {
                for (int d : { X, Y, Z }) {
                    const int comp = (ft == E ? 0 : 3) + d;
                    out[static_cast<std::size_t>(comp * ndofs + i)] +=
                        scale * evalDirectPlanewaveComponent(src, x, y, z, time, ft, d);
                }
            }
        }
    }

#ifdef SEMBA_DGTD_ENABLE_CUDA
    // HostRead above can leave these vectors host-authoritative; Mult's GPU
    // kernel needs device-valid pointers on the next eval_tfsf_planewave_gpu.
    if (mfem::Device::Allows(mfem::Backend::CUDA)) {
        (void)const_cast<mfem::Vector&>(dof_coords_).Read();
        (void)const_cast<mfem::Vector&>(tfsf_sign_).Read();
        (void)const_cast<mfem::Vector&>(direct_source_params_).Read();
        (void)const_cast<mfem::Array<int>&>(direct_source_meta_).Read();
    }
#endif
}

int SourcesManager::findTFSFDofAtPosition(const Source::Position& pos, double tol) const
{
    if (!direct_eval_ready_ || dof_coords_.Size() == 0) {
        return -1;
    }
    const int dim = pos.Size();
    const int ndofs = dof_coords_.Size() / 3;
    const double* coords = dof_coords_.HostRead();
    for (int i = 0; i < ndofs; ++i) {
        double dist2 = 0.0;
        for (int d = 0; d < dim && d < 3; ++d) {
            const double diff = pos[d] - coords[3 * i + d];
            dist2 += diff * diff;
        }
        if (dist2 <= tol * tol) {
            return i;
        }
    }
    return -1;
}

}