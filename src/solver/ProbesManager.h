#pragma once

#include <iostream>
#include <fstream>

#include "components/Probes.h"
#include "components/SubMesher.h"
#include "components/RCSSurfaceExporter.h"
#include "solver/SolverOptions.h"

namespace maxwell {

class SourcesManager;  // Forward declaration

std::string getRunModeTag();

class NearFieldReqs {
public:

    NearFieldReqs(const NearFieldProbe&, const mfem::DG_FECollection* fec, mfem::ParFiniteElementSpace& fes, Fields<ParFiniteElementSpace, ParGridFunction>&);

    mfem::SubMesh* getSubMesh() { return ntff_smsh_.getSubMesh(); }
    const mfem::GridFunction& getConstField(const FieldType& f, const Direction& d) const { return fields_.get(f, d); }
    mfem::GridFunction& getConstField(const FieldType& f, const Direction& d) { return fields_.get(f, d); }
    void updateFields();

private:

    NearToFarFieldSubMesher ntff_smsh_;
    std::unique_ptr<mfem::FiniteElementSpace> sfes_;
    Fields<FiniteElementSpace, GridFunction> fields_;
    Fields<ParFiniteElementSpace, ParGridFunction>& gFields_;
    TransferMaps tMaps_;

};

class ProbesManager {
public:
    ProbesManager() = delete;
    ProbesManager(Probes, mfem::ParFiniteElementSpace&, Fields<ParFiniteElementSpace, ParGridFunction>&, const SolverOptions&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(Time);
    /// True when at least one probe will read host field data this cycle.
    bool needsHostSyncThisStep(Time) const;
    void recalculateExportSteps(double dt);
    void setFinalTime(double final_time);

    const FieldProbe& getFieldProbe(const std::size_t i) const;
    const PointProbe& getPointProbe(const std::size_t i) const;

    void setCaseName(const std::string name) {
        caseName_ = name;
        initRCSSurfaceExporters();
    }
    void initPointFieldProbeExport();
    void initTFSFExport(SourcesManager* srcmngr, const mfem::Array<int>* tfsf_mapping) {
        srcmngr_ = srcmngr;
        tfsf_mapping_ = tfsf_mapping;
    }
    void printTimingSummaryAndReset() const;

    Probes probes;

private:
    struct TimingStats {
        double exporter_ms{0.0};
        double field_ms{0.0};
        double point_ms{0.0};
        double nearfield_ms{0.0};
        double snapshot_ms{0.0};
        double rcs_ms{0.0};
        double mor_ms{0.0};
        int update_calls{0};
    };

    struct FESPoint {
        int elementId;
        mfem::IntegrationPoint iP;
    };

    struct PointProbeCollection {
        FESPoint fesPoint;
        const mfem::GridFunction& field_Ex;
        const mfem::GridFunction& field_Ey;
        const mfem::GridFunction& field_Ez;
        const mfem::GridFunction& field_Hx;
        const mfem::GridFunction& field_Hy;
        const mfem::GridFunction& field_Hz;
    };

    struct FieldProbeCollection {
        FESPoint fesPoint;
        const mfem::GridFunction& field;
    };

    struct ExporterContext {
        int save_count{0};
        double next_save_time{0.0};
        double dt_save{0.0};
        bool initialized{false};
    };

    int cycle_{ 0 };
    double finalTime_;

    std::map<const ExporterProbe*, ExporterContext> exporterContexts_;
    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointProbe*, PointProbeCollection> pointProbesCollection_;
    std::map<const FieldProbe*, FieldProbeCollection> fieldProbesCollection_;
    std::map<const NearFieldProbe*, DataCollection> nearFieldProbesCollection_;
    std::map<const DomainSnapshotProbe*, DomainSnapshotDataCollection> domainSnapshotProbesCollection_;

    std::string caseName_;
    
    mfem::ParFiniteElementSpace& fes_;
    Fields<ParFiniteElementSpace, ParGridFunction>* fields_;

    std::map<const NearFieldProbe*, std::unique_ptr<NearFieldReqs>> nearFieldReqs_;
    std::map<const RCSSurfaceProbe*, std::unique_ptr<RCSSurfaceExporter>> rcsSurfaceExporters_;
    std::map<int, std::ofstream> fieldProbeFiles_;
    std::map<int, std::ofstream> pointProbeFiles_;

    struct MORStateContext {
        int save_count{0};
        double next_save_time{0.0};
        double dt_save{0.0};
        bool initialized{false};
        std::string export_dir;
    };
    std::map<const MORStateProbe*, MORStateContext> morStateContexts_;
    
    SourcesManager* srcmngr_{nullptr};
    const mfem::Array<int>* tfsf_mapping_{nullptr};
    bool is_sgbc_solver_{false};
    mutable TimingStats timingStats_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields<ParFiniteElementSpace, ParGridFunction>&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields<ParFiniteElementSpace, ParGridFunction>&) const;
    FieldProbeCollection buildFieldProbeCollectionInfo(const FieldProbe&, Fields<ParFiniteElementSpace, ParGridFunction>&) const;
    DataCollection buildNearFieldDataCollectionInfo(const NearFieldProbe&, Fields<ParFiniteElementSpace, ParGridFunction>&) const;
    DomainSnapshotDataCollection buildDomainSnapshotDataCollection(const DomainSnapshotProbe& p, Fields<ParFiniteElementSpace, ParGridFunction>& fields) const;

    void updateProbe(ExporterProbe&, Time);
    void updateProbe(FieldProbe&, Time);
    void updateProbe(PointProbe&, Time);
    void updateProbe(NearFieldProbe&, Time);
    void updateProbe(DomainSnapshotProbe&, Time);
    void updateProbe(RCSSurfaceProbe&, Time);
    void updateProbe(MORStateProbe&, Time);
    void initRCSSurfaceExporters();
};


}