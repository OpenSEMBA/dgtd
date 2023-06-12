#pragma once

#include "components/Probes.h"
#include "evolution/Fields.h"
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

    void updateProbes(Time);

    const PointProbe& getPointProbe(const std::size_t i) const;
    const FieldProbe& getFieldProbe(const std::size_t i) const;
    //const EnergyProbe& getEnergyProbe(const std::size_t i) const;

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
        const Fields& fields;
    };

    //struct EnergyProbeCollection {
    //    FiniteElementSpace fes;
    //    const Fields& fields;
    //};

    int cycle_{ 0 };

    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointProbe*, PointProbeCollection> pointProbesCollection_;
    std::map<const FieldProbe*, FieldProbeCollection> fieldProbesCollection_;
    //std::map<const EnergyProbe*, EnergyProbeCollection> energyProbesCollection_;
    
    const mfem::FiniteElementSpace& fes_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields&) const;
    FieldProbeCollection buildFieldProbeCollectionInfo(const FieldProbe&, Fields&) const;
    //EnergyProbeCollection buildEnergyProbeCollectionInfo(const mfem::FiniteElementSpace& fes, const Fields&) const;

    void updateProbe(ExporterProbe&, Time);
    void updateProbe(PointProbe&, Time);
    void updateProbe(FieldProbe&, Time);
    //void updateProbe(EnergyProbe&, Time);
};

}