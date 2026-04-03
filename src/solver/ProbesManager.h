#pragma once

#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "general/text.hpp"

#include "components/Probes.h"
#include "components/SubMesher.h"
#include "components/RCSSurfaceExporter.h"
#include "solver/SolverOptions.h"

namespace maxwell {

Array<int> buildSurfaceMarker(const std::vector<int>& tags, const ParFiniteElementSpace&);
std::string getRunModeTag();

struct TransferMaps {

    mfem::TransferMap tMapEx;
    mfem::TransferMap tMapEy;
    mfem::TransferMap tMapEz;
    mfem::TransferMap tMapHx;
    mfem::TransferMap tMapHy;
    mfem::TransferMap tMapHz;

    TransferMaps(Fields<ParFiniteElementSpace, ParGridFunction>& src, Fields<FiniteElementSpace, GridFunction>& dst) :
        tMapEx{ mfem::TransferMap(src.get(E, X), dst.get(E, X)) },
        tMapEy{ mfem::TransferMap(src.get(E, Y), dst.get(E, Y)) },
        tMapEz{ mfem::TransferMap(src.get(E, Z), dst.get(E, Z)) },
        tMapHx{ mfem::TransferMap(src.get(H, X), dst.get(H, X)) },
        tMapHy{ mfem::TransferMap(src.get(H, Y), dst.get(H, Y)) },
        tMapHz{ mfem::TransferMap(src.get(H, Z), dst.get(H, Z)) }
    {}

    void transferFields(const Fields<ParFiniteElementSpace, ParGridFunction>& src, Fields<FiniteElementSpace, GridFunction>& dst)
    {
        tMapEx.Transfer(src.get(E, X), dst.get(E, X));
        tMapEy.Transfer(src.get(E, Y), dst.get(E, Y));
        tMapEz.Transfer(src.get(E, Z), dst.get(E, Z));
        tMapHx.Transfer(src.get(H, X), dst.get(H, X));
        tMapHy.Transfer(src.get(H, Y), dst.get(H, Y));
        tMapHz.Transfer(src.get(H, Z), dst.get(H, Z));
    }
};

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
    void recalculateExportSteps(double dt);

    const FieldProbe& getFieldProbe(const std::size_t i) const;
    const PointProbe& getPointProbe(const std::size_t i) const;

    void setCaseName(const std::string name) {
        caseName_ = name;
        initRCSSurfaceExporters();
    }
    void initPointFieldProbeExport();

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
    std::map<const DomainSnapshotProbe*, DomainSnapshotDataCollection> domainSnapshotProbesCollection_;

    std::string caseName_;
    
    mfem::ParFiniteElementSpace& fes_;
    Fields<ParFiniteElementSpace, ParGridFunction>* fields_;

    std::map<const NearFieldProbe*, std::unique_ptr<NearFieldReqs>> nearFieldReqs_;
    std::map<const RCSSurfaceProbe*, std::unique_ptr<RCSSurfaceExporter>> rcsSurfaceExporters_;
    std::map<int, std::ofstream> fieldProbeFiles_;
    std::map<int, std::ofstream> pointProbeFiles_;
    
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
    void initRCSSurfaceExporters();
};


}