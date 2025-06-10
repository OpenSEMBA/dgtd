#pragma once

#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "general/text.hpp"

#include "components/Probes.h"
#include "components/SubMesher.h"
#include "evolution/Fields.h"
#include "solver/SolverOptions.h"

namespace maxwell {

Array<int> buildSurfaceMarker(const std::vector<int>& tags, const ParFiniteElementSpace&);

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

class NearFieldReqs {
public:

    NearFieldReqs(const NearFieldProbe&, const mfem::DG_FECollection* fec, mfem::ParFiniteElementSpace& fes, Fields&);

    mfem::ParSubMesh* getSubMesh() { return ntff_smsh_.getSubMesh(); }
    const mfem::ParGridFunction& getConstField(const FieldType& f, const Direction& d) { return fields_.get(f, d); }
    mfem::ParGridFunction& getField(const FieldType& f, const Direction& d) { return fields_.get(f, d); }
    void updateFields();

private:

    void assignGlobalFieldsReferences(Fields& global);

    NearToFarFieldSubMesher ntff_smsh_;
    std::unique_ptr<mfem::ParFiniteElementSpace> sfes_;
    Fields fields_;
    Fields& gFields_;
    TransferMaps tMaps_;

};

class ProbesManager {
public:
    ProbesManager() = delete;
    ProbesManager(Probes, mfem::ParFiniteElementSpace&, Fields&, const SolverOptions&);
    
    ProbesManager(const ProbesManager&) = delete;
    ProbesManager(ProbesManager&&) = default;
    ~ProbesManager() = default;
    ProbesManager& operator=(const ProbesManager&) = delete;
    ProbesManager& operator=(ProbesManager&&) = default;

    void updateProbes(Time);

    const FieldProbe& getFieldProbe(const std::size_t i) const;
    const PointProbe& getPointProbe(const std::size_t i) const;

    Probes probes;

private:

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

    int cycle_{ 0 };
    double finalTime_;

    std::map<const ExporterProbe*, mfem::ParaViewDataCollection> exporterProbesCollection_;
    std::map<const PointProbe*, PointProbeCollection> pointProbesCollection_;
    std::map<const FieldProbe*, FieldProbeCollection> fieldProbesCollection_;
    std::map<const NearFieldProbe*, DataCollection> nearFieldProbesCollection_;
    
    mfem::ParFiniteElementSpace& fes_;
    Fields* fields_;

    std::map<const NearFieldProbe*, std::unique_ptr<NearFieldReqs>> nearFieldReqs_;
    
    mfem::ParaViewDataCollection buildParaviewDataCollectionInfo(const ExporterProbe&, Fields&) const;
    PointProbeCollection buildPointProbeCollectionInfo(const PointProbe&, Fields&) const;
    FieldProbeCollection buildFieldProbeCollectionInfo(const FieldProbe&, Fields&) const;
    DataCollection buildNearFieldDataCollectionInfo(const NearFieldProbe&, Fields&) const;

    void updateProbe(ExporterProbe&, Time);
    void updateProbe(FieldProbe&, Time);
    void updateProbe(PointProbe&, Time);
    void updateProbe(NearFieldProbe&, Time);
};


}