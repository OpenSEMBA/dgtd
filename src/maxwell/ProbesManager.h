#pragma once

#include "Probes.h"

namespace maxwell {

class ProbesManager {
public:
    ProbesManager() = default;
    ProbesManager(Probes, const mfem::FiniteElementSpace*, FieldViews&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(double time);

    const PointsProbe* getPointsProbe(const std::size_t i) const;

private:
    struct FESPoint {
        int elementId;
        mfem::IntegrationPoint iP;
    };

    struct PointsProbeCollection {
        std::vector<FESPoint> fesPoints;
        const mfem::GridFunction* field;
    };

    int cycle_{ 0 };

    Probes probes_;
    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointsProbe*, PointsProbeCollection> pointProbesCollection_;
    
    const mfem::FiniteElementSpace* fes_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollection(FieldViews& fields) const;
    PointsProbeCollection buildPointsProbeCollection(const PointsProbe&, FieldViews& fields) const;
    
    void updateProbe(ExporterProbe&, double time);
    void updateProbe(PointsProbe&, double time);
};

}