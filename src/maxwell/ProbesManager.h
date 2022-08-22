#pragma once

#include "Probes.h"

namespace maxwell {

class ProbesManager {
public:
    ProbesManager() = default;
    ProbesManager(Probes, const mfem::FiniteElementSpace*, const FieldViews&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(bool done, int cycle);

    const PointsProbe* getPointsProbe(const std::size_t i) const;
    
    using IntegrationPoint = mfem::IntegrationPoint;
    using IntegrationPointsSet = std::vector<std::vector<IntegrationPoint>>;

    struct PointsProbeInFES {
        mfem::Array<int> elemIds;
        IntegrationPointsSet integPointSet;
        const GridFunction& field;
    };

private:

    int vis_steps = 1;
    Probes probes_;
    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> dataCollection_;
    std::map<const PointsProbe*, PointsProbeInFES> probesToFES_;
    
    PointsProbeInFES buildElemAndIntegrationPointArrays(const PointsProbe&) const;
    FieldFrame getFieldForPointsProbe(const PointsProbe& p) const;

    void storeInitialVisualizationValues();

};

}