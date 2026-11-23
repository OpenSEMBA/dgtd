#pragma once

#include "Model.h"

#include <mfem.hpp>
#include "mesh/submesh/psubmesh.hpp"

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
                          const std::vector<std::string>& requested,
                          bool smooth_plasma = false);

    bool empty() const { return fields_.empty(); }
    bool isTimeDependent() const;
    const std::vector<Field>& fields() const { return fields_; }
    const std::vector<Field>& plasmaFields() const { return plasma_fields_; }
    const std::vector<Field>& plasmaCellFields() const { return plasma_cell_fields_; }
    mfem::ParMesh* plasmaMesh() { return plasma_submesh_.get(); }
    int plasmaOrder() const { return plasma_fec_ ? plasma_fec_->GetOrder() : 0; }
    bool hasPlasmaFields() const { return !plasma_fields_.empty(); }
    bool hasPlasmaCellFields() const { return !plasma_cell_fields_.empty(); }

    void registerFields(mfem::ParaViewDataCollection& pd) const;
    void registerPlasmaFields(mfem::ParaViewDataCollection& pd) const;
    void registerPlasmaCellFields(mfem::ParaViewDataCollection& pd) const;

    static const std::vector<std::string>& allFieldNames();
    static std::vector<std::string> resolveRequested(
        const std::vector<std::string>& requested);
    /// Fields that make sense for this model (no PML fields without PML,
    /// no Debye fields without Debye, etc.).
    static std::vector<std::string> applicableFieldNames(const Model& model);

private:
    const Model& model_;
    mfem::ParFiniteElementSpace& fes_;
    bool smooth_plasma_{false};
    std::vector<Field> fields_;
    std::unique_ptr<mfem::ParSubMesh> plasma_submesh_;
    std::unique_ptr<mfem::H1_FECollection> plasma_fec_;
    std::unique_ptr<mfem::ParFiniteElementSpace> plasma_fes_;
    std::vector<Field> plasma_fields_;
    std::unique_ptr<mfem::L2_FECollection> plasma_cell_fec_;
    std::unique_ptr<mfem::ParFiniteElementSpace> plasma_cell_fes_;
    std::vector<Field> plasma_cell_fields_;

    void build();
    void buildPlasmaVisualization();
    mfem::ParGridFunction* findField(const std::string& name);
};

} // namespace maxwell
