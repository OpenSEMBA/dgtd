#pragma once

#include <iostream>
#include <fstream>
#include <mfem.hpp>

#ifdef __linux__
#     include <mfem/general/text.hpp>
#elif _WIN32
#     include <general/text.hpp>
#endif

#include <components/Probes.h>
#include <evolution/Fields.h>
#include <solver/SolverOptions.h>
#include <components/SubMesher.h>

namespace maxwell {

struct TransferMaps {

    mfem::TransferMap tMapEx;
    mfem::TransferMap tMapEy;
    mfem::TransferMap tMapEz;
    mfem::TransferMap tMapHx;
    mfem::TransferMap tMapHy;
    mfem::TransferMap tMapHz;

    TransferMaps(Fields& src, Fields& dst) :
        tMapEx{ mfem::TransferMap(src.get(E, X), dst.get(E, X)) },
        tMapEy{ mfem::TransferMap(src.get(E, Y), dst.get(E, Y)) },
        tMapEz{ mfem::TransferMap(src.get(E, Z), dst.get(E, Z)) },
        tMapHx{ mfem::TransferMap(src.get(H, X), dst.get(H, X)) },
        tMapHy{ mfem::TransferMap(src.get(H, Y), dst.get(H, Y)) },
        tMapHz{ mfem::TransferMap(src.get(H, Z), dst.get(H, Z)) }
    {}

    void transferFields(const Fields&, Fields&);
};

class NearToFarFieldReqs {
public:

    NearToFarFieldReqs(const NearToFarFieldProbe&, const mfem::DG_FECollection* fec, mfem::FiniteElementSpace& fes, Fields&);

    mfem::SubMesh* getSubMesh() { return ntff_smsh_.getSubMesh(); }
    const mfem::GridFunction& getConstField(const FieldType& f, const Direction& d) { return fields_.get(f, d); }
    mfem::GridFunction& getField(const FieldType& f, const Direction& d) { return fields_.get(f, d); }
    void updateFields();

private:

    void assignGlobalFieldsReferences(Fields& global);

    NearToFarFieldSubMesher ntff_smsh_;
    std::unique_ptr<mfem::FiniteElementSpace> sfes_;
    Fields fields_;
    Fields& gFields_;
    TransferMaps tMaps_;

};

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
    std::map<const NearToFarFieldProbe*, DataCollection> nearToFarFieldProbesCollection_;
    
    mfem::FiniteElementSpace& fes_;

    std::map<const NearToFarFieldProbe*, std::unique_ptr<NearToFarFieldReqs>> nearToFarFieldReqs_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields&) const;
    FieldProbeCollection buildFieldProbeCollectionInfo(const FieldProbe&, Fields&) const;
    DataCollection buildNearToFarFieldDataCollectionInfo(const NearToFarFieldProbe&, Fields&) const;

    void updateProbe(ExporterProbe&, Time);
    void updateProbe(PointProbe&, Time);
    void updateProbe(FieldProbe&, Time);
    void updateProbe(NearToFarFieldProbe&, Time);
};


}