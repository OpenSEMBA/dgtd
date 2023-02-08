#pragma once

#include "Probes.h"
#include "Fields.h"
#include "SolverOptions.h"

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

    void updateProbes(double time);

    const PointProbe& getPointProbe(const std::size_t i) const;

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

    int cycle_{ 0 };

    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointProbe*, PointProbeCollection> pointProbesCollection_;
    
    const mfem::FiniteElementSpace& fes_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields&) const;

    void updateProbe(ExporterProbe&, double time);
    void updateProbe(PointProbe&, double time);
};

}