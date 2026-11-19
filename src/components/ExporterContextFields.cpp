#include "ExporterContextFields.h"

#include "math/PhysicalConstants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>

namespace maxwell {

namespace {

const std::vector<std::string> kAllFields = {
    "eps_r", "mu_r", "sigma_S_per_m", "mat_tag", "h_m",
    "eps_inf", "fp_Hz", "f1_Hz", "nuc_Hz", "eps_s", "tau_s",
    "pml_sigma_x", "pml_sigma_y", "pml_sigma_z",
    "pml_kappa_x", "pml_kappa_y", "pml_kappa_z",
    "pml_alpha_x", "pml_alpha_y", "pml_alpha_z",
    "pml_depth_x", "pml_depth_y", "pml_depth_z",
    "pml_mask"
};

double defaultValueFor(const std::string& name)
{
    if (name == "eps_r" || name == "mu_r" || name == "eps_inf" ||
        name == "eps_s" || name.rfind("pml_kappa_", 0) == 0) {
        return 1.0;
    }
    return 0.0;
}

} // namespace

const std::vector<std::string>& ExporterContextFields::allFieldNames()
{
    return kAllFields;
}

std::vector<std::string> ExporterContextFields::resolveRequested(
    const std::vector<std::string>& requested)
{
    if (requested.empty()) {
        return kAllFields;
    }
    if (requested.size() == 1) {
        if (requested[0] == "all") {
            return kAllFields;
        }
        if (requested[0] == "none") {
            return {};
        }
    }

    std::vector<std::string> resolved;
    resolved.reserve(requested.size());
    for (const auto& name : requested) {
        if (std::find(kAllFields.begin(), kAllFields.end(), name) ==
            kAllFields.end()) {
            throw std::runtime_error(
                "Unknown exporter context field '" + name + "'. Expected one of: "
                + [] {
                      std::string joined;
                      for (const auto& f : kAllFields) {
                          if (!joined.empty()) joined += ", ";
                          joined += f;
                      }
                      return joined;
                  }());
        }
        resolved.push_back(name);
    }
    return resolved;
}

std::vector<std::string> ExporterContextFields::applicableFieldNames(
    const Model& model)
{
    std::vector<std::string> names = {"eps_r", "mu_r"};

    bool any_sigma = false;
    for (const auto& entry : model.getGeomTagToMaterial()) {
        if (entry.second.getConductivity() != 0.0) {
            any_sigma = true;
            break;
        }
    }
    if (any_sigma) {
        names.push_back("sigma_S_per_m");
    }

    names.push_back("mat_tag");
    names.push_back("h_m");

    if (model.hasLorentz()) {
        names.push_back("eps_inf");
        names.push_back("fp_Hz");
        bool has_resonance = false;
        for (const auto& pole : model.getLorentzProperties()) {
            if (pole.omega_1_profile != nullptr || pole.omega_1 != 0.0) {
                has_resonance = true;
                break;
            }
        }
        if (has_resonance) {
            names.push_back("f1_Hz");
        }
        names.push_back("nuc_Hz");
    }

    if (model.hasDebye()) {
        if (!model.hasLorentz()) {
            names.push_back("eps_inf");
        }
        names.push_back("eps_s");
        names.push_back("tau_s");
    }

    if (model.hasPML()) {
        std::set<Direction> active;
        for (const auto& pml : model.getPMLProperties()) {
            active.insert(pml.active_axes.begin(), pml.active_axes.end());
        }
        auto add_axis_fields = [&](const std::string& prefix) {
            for (Direction axis : active) {
                const char* suffix =
                    axis == X ? "x" : (axis == Y ? "y" : "z");
                names.push_back(prefix + suffix);
            }
        };
        add_axis_fields("pml_sigma_");
        add_axis_fields("pml_kappa_");
        add_axis_fields("pml_alpha_");
        add_axis_fields("pml_depth_");
        names.push_back("pml_mask");
    }

    return names;
}

ExporterContextFields::ExporterContextFields(
    const Model& model,
    mfem::ParFiniteElementSpace& fes,
    const std::vector<std::string>& requested)
    : model_(model), fes_(fes)
{
    auto names = resolveRequested(requested);
    if (requested.empty() ||
        (requested.size() == 1 && requested[0] == "all")) {
        names = applicableFieldNames(model);
    }

    for (const auto& name : names) {
        Field field;
        field.name = name;
        field.gf = std::make_unique<mfem::ParGridFunction>(&fes_);
        *field.gf = defaultValueFor(name);
        // On device builds the scalar assignment above may only initialize the
        // device copy.  Pull the defaults to the host before build() edits them.
        field.gf->HostRead();
        fields_.push_back(std::move(field));
    }
    build();
}

bool ExporterContextFields::isTimeDependent() const
{
    return std::any_of(fields_.begin(), fields_.end(),
                       [](const Field& f) { return f.time_dependent; });
}

void ExporterContextFields::registerFields(mfem::ParaViewDataCollection& pd) const
{
    for (const auto& field : fields_) {
        pd.RegisterField(field.name, field.gf.get());
    }
}

mfem::ParGridFunction* ExporterContextFields::findField(const std::string& name)
{
    for (auto& field : fields_) {
        if (field.name == name) {
            return field.gf.get();
        }
    }
    return nullptr;
}

void ExporterContextFields::build()
{
    const auto& mesh = *fes_.GetMesh();
    const double c_si = physicalConstants::speedOfLight_SI;
    const double two_pi = 2.0 * M_PI;
    const double z0_si = physicalConstants::freeSpaceImpedance_SI;

    auto* eps_r = findField("eps_r");
    auto* mu_r = findField("mu_r");
    auto* sigma = findField("sigma_S_per_m");
    auto* mat_tag = findField("mat_tag");
    auto* h_m = findField("h_m");
    auto* eps_inf = findField("eps_inf");
    auto* fp = findField("fp_Hz");
    auto* f1 = findField("f1_Hz");
    auto* nuc = findField("nuc_Hz");
    auto* eps_s = findField("eps_s");
    auto* tau = findField("tau_s");
    auto* pml_mask = findField("pml_mask");

    const std::array<mfem::ParGridFunction*, 3> pml_sigma{{
        findField("pml_sigma_x"), findField("pml_sigma_y"),
        findField("pml_sigma_z")}};
    const std::array<mfem::ParGridFunction*, 3> pml_kappa{{
        findField("pml_kappa_x"), findField("pml_kappa_y"),
        findField("pml_kappa_z")}};
    const std::array<mfem::ParGridFunction*, 3> pml_alpha{{
        findField("pml_alpha_x"), findField("pml_alpha_y"),
        findField("pml_alpha_z")}};
    const std::array<mfem::ParGridFunction*, 3> pml_depth{{
        findField("pml_depth_x"), findField("pml_depth_y"),
        findField("pml_depth_z")}};

    const auto& materials = model_.getGeomTagToMaterial();
    const auto* pml_profiles = model_.getPMLProfileData();

    mfem::Array<int> dofs;
    mfem::Vector x(mesh.SpaceDimension());

    for (int el = 0; el < mesh.GetNE(); ++el) {
        const int attr = mesh.GetAttribute(el);
        fes_.GetElementDofs(el, dofs);

        const auto* fe = fes_.GetFE(el);
        auto* tr = fes_.GetElementTransformation(el);
        if (fe->GetNodes().Size() != dofs.Size()) {
            throw std::runtime_error(
                "ExporterContextFields: element " + std::to_string(el)
                + " has " + std::to_string(dofs.Size()) + " DOFs but "
                + std::to_string(fe->GetNodes().Size())
                + " nodal points; cannot attach context fields per DOF.");
        }

        const Material* mat = nullptr;
        const auto mat_it = materials.find(attr);
        if (mat_it != materials.end()) {
            mat = &mat_it->second;
        }
        const auto* lorentz = model_.findLorentz(attr);
        const auto* debye = model_.findDebye(attr);
        const bool is_pml = model_.isPMLAttribute(attr);
        const double h = fes_.GetMesh()->GetElementSize(el);

        for (int j = 0; j < dofs.Size(); ++j) {
            const int dof = dofs[j] < 0 ? -1 - dofs[j] : dofs[j];
            const auto& ip = fe->GetNodes().IntPoint(j);
            tr->Transform(ip, x);

            if (mat_tag) { (*mat_tag)(dof) = attr; }
            if (h_m) { (*h_m)(dof) = h; }

            if (mat) {
                if (eps_r) { (*eps_r)(dof) = mat->getPermittivity(); }
                if (mu_r) { (*mu_r)(dof) = mat->getPermeability(); }
                if (sigma) {
                    (*sigma)(dof) = mat->getConductivity() / z0_si;
                }
            }

            if (lorentz) {
                if (eps_inf) { (*eps_inf)(dof) = lorentz->eps_inf; }
                if (fp) {
                    (*fp)(dof) = lorentz->omega_p_gf
                        ? (*lorentz->omega_p_gf)(dof) * c_si / two_pi
                        : lorentz->evalOmegaP(x) * c_si / two_pi;
                }
                if (f1) {
                    (*f1)(dof) = lorentz->omega_1_gf
                        ? (*lorentz->omega_1_gf)(dof) * c_si / two_pi
                        : lorentz->evalOmega1(x) * c_si / two_pi;
                }
                if (nuc) {
                    (*nuc)(dof) = lorentz->gamma_gf
                        ? (*lorentz->gamma_gf)(dof) * c_si / two_pi
                        : lorentz->evalGamma(x) * c_si / two_pi;
                }
            }

            if (debye) {
                if (eps_inf) { (*eps_inf)(dof) = debye->eps_inf; }
                if (eps_s) { (*eps_s)(dof) = debye->eps_s; }
                if (tau) { (*tau)(dof) = debye->tau_solver / c_si; }
            }

            if (is_pml && pml_profiles != nullptr) {
                if (pml_mask) { (*pml_mask)(dof) = 1.0; }
                for (int d = 0; d < 3; ++d) {
                    PMLDirectionProfiles profile;
                    pml_profiles->evaluateAtTransform(
                        *tr, ip, static_cast<Direction>(d), profile);
                    if (pml_sigma[d]) { (*pml_sigma[d])(dof) = profile.sigma; }
                    if (pml_kappa[d]) { (*pml_kappa[d])(dof) = profile.kappa; }
                    if (pml_alpha[d]) { (*pml_alpha[d])(dof) = profile.alpha; }
                    if (pml_depth[d]) { (*pml_depth[d])(dof) = profile.depth; }
                }
            }
        }
    }

    // The values above are written through the host pointer without touching
    // MFEM's memory validity flags.  Mark the host side as authoritative so a
    // CUDA build uploads these fields instead of exporting the initial defaults.
    for (auto& field : fields_) {
        field.gf->HostWrite();
    }
}

} // namespace maxwell
