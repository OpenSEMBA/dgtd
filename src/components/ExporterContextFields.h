#pragma once

#include "Model.h"

#include <mfem.hpp>

#include <memory>
#include <string>
#include <vector>

namespace maxwell {

/// Builds the static, DOF-attached context fields written by the ParaView
/// exporter (material constants, dispersive plasma parameters, PML profiles).
///
/// Every field is a scalar ParGridFunction on the same DG space as the E/H
/// fields. Values are evaluated at each DOF's physical coordinate, so spatial
/// profiles (e.g. omega_p / nu_c gridfunctions) vary DOF by DOF rather than
/// being element-constant.
class ExporterContextFields {
public:
    struct Field {
        std::string name;
        std::unique_ptr<mfem::ParGridFunction> gf;
        /// True when the value can change during the run (none today; reserved
        /// for future fields such as P/J auxiliaries).
        bool time_dependent{false};
    };

    ExporterContextFields(const Model& model,
                          mfem::ParFiniteElementSpace& fes,
                          const std::vector<std::string>& requested);

    bool empty() const { return fields_.empty(); }
    bool isTimeDependent() const;
    const std::vector<Field>& fields() const { return fields_; }

    void registerFields(mfem::ParaViewDataCollection& pd) const;

    static const std::vector<std::string>& allFieldNames();
    static std::vector<std::string> resolveRequested(
        const std::vector<std::string>& requested);
    /// Fields that make sense for this model (no PML fields without PML,
    /// no Debye fields without Debye, etc.).
    static std::vector<std::string> applicableFieldNames(const Model& model);

private:
    const Model& model_;
    mfem::ParFiniteElementSpace& fes_;
    std::vector<Field> fields_;

    void build();
    mfem::ParGridFunction* findField(const std::string& name);
};

} // namespace maxwell
