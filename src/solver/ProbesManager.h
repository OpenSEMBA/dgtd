#pragma once

#include <iostream>
#include <fstream>
#include <mfem.hpp>

#include <components/Probes.h>
#include <evolution/Fields.h>
#include <solver/SolverOptions.h>
#include <components/Exporter.h>
#include <components/SubMesher.h>

namespace maxwell {

class ProbesManager {
public:
    ProbesManager() = delete;
    ProbesManager(Probes, mfem::FiniteElementSpace&, Fields&, const SolverOptions&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(Time);

    const PointProbe& getPointProbe(const std::size_t i) const;
    const FieldProbe& getFieldProbe(const std::size_t i) const;

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

    int cycle_{ 0 };
    double finalTime_;

    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointProbe*, PointProbeCollection> pointProbesCollection_;
    std::map<const FieldProbe*, FieldProbeCollection> fieldProbesCollection_;
    std::map<const NearToFarFieldProbe*, NearToFarFieldDataCollection> nearToFarFieldProbesCollection_;
    
    mfem::FiniteElementSpace& fes_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields&) const;
    FieldProbeCollection buildFieldProbeCollectionInfo(const FieldProbe&, Fields&) const;
    NearToFarFieldDataCollection buildNearToFarFieldDataCollectionInfo(const NearToFarFieldProbe&, const DG_FECollection&, Fields&) const;

    void performNearToFarFieldExports(const NearToFarFieldProbe&, NearToFarFieldSubMesher&);

    void updateProbe(ExporterProbe&, Time);
    void updateProbe(PointProbe&, Time);
    void updateProbe(FieldProbe&, Time);
    void updateProbe(NearToFarFieldProbe&, Time);
};

}