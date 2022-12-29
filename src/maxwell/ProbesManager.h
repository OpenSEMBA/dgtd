#pragma once

#include "Probes.h"
#include "Fields.h"

namespace maxwell {

class ProbesManager {
public:
    ProbesManager() = delete;
    ProbesManager(Probes, const mfem::FiniteElementSpace&, Fields&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(double time);

    const PointsProbe& getPointsProbe(const std::size_t i) const;

private:
    struct FESPoint {
        int elementId;
        mfem::IntegrationPoint iP;
    };

    struct PointsProbeCollection {
        std::vector<FESPoint> fesPoints;
        const mfem::GridFunction& field;
    };

    int cycle_{ 0 };

    Probes probes_;
    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointsProbe*, PointsProbeCollection> pointProbesCollection_;
    
    const mfem::FiniteElementSpace& fes_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollection(const ExporterProbe&, Fields&) const;
    PointsProbeCollection buildPointsProbeCollection(const PointsProbe&, Fields&) const;
    
    void updateProbe(ExporterProbe&, double time);
    void updateProbe(PointsProbe&, double time);
};

}