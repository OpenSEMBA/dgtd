#pragma once

#include "components/Probes.h"
#include "evolution/Fields.h"
#include "SolverOptions.h"
#include <components/Exporter.h>

namespace maxwell {

class ProbesManager {
public:
    ProbesManager() = delete;
    ProbesManager(Probes, const mfem::FiniteElementSpace&, Fields&, const SolverOptions&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(Time);

    const PointProbe& getPointProbe(const std::size_t i) const;
    const FieldProbe& getFieldProbe(const std::size_t i) const;

    void initNearToFarFieldProbeDataCollection(NearToFarFieldProbe&, Fields&);

    Probes probes;

private:

    struct FESPoint {
        int elementId;
        mfem::IntegrationPoint iP;
    };

    struct PointProbeCollection {
        FESPoint fesPoint;
        const mfem::GridFunction& field;
    };

    struct FieldProbeCollection {
        FESPoint fesPoint;
        const mfem::GridFunction& field_Ex;
        const mfem::GridFunction& field_Ey;
        const mfem::GridFunction& field_Ez;
        const mfem::GridFunction& field_Hx;
        const mfem::GridFunction& field_Hy;
        const mfem::GridFunction& field_Hz;
    };

    using NearToFarFieldProbeCollection = NearToFarFieldExporter;

    int cycle_{ 0 };
    double finalTime_;

    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointProbe*, PointProbeCollection> pointProbesCollection_;
    std::map<const FieldProbe*, FieldProbeCollection> fieldProbesCollection_;
    std::map<const NearToFarFieldProbe*, NearToFarFieldExporter> nearToFarFieldProbesCollection_;
    
    const mfem::FiniteElementSpace& fes_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields&) const;
    FieldProbeCollection buildFieldProbeCollectionInfo(const FieldProbe&, Fields&) const;
    NearToFarFieldExporter buildNearToFarFieldDataCollectionInfo(const NearToFarFieldProbe&, Fields&);


    void updateProbe(ExporterProbe&, Time);
    void updateProbe(PointProbe&, Time);
    void updateProbe(FieldProbe&, Time);
    void updateProbe(NearToFarFieldProbe&, Time);
};

}