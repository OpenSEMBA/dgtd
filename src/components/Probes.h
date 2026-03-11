#pragma once

#include <mfem.hpp>

#include "Types.h"
#include "evolution/Fields.h"

namespace maxwell {

struct ExporterProbe {
    std::string name{"MaxwellView"};
    int visSteps{ 10 };
};

struct NearFieldProbe {
    std::string name = std::string("NearFieldProbe");
    std::string exportPath = std::string("./defaultExportPath/" + name + "/");
    int expSteps{ 10 };
    std::vector<int> tags;
};

struct DomainSnapshotProbe {
    std::string name = std::string("DomainSnapshot");
    int expSteps { 10 };
};

struct DomainSnapshotDataCollection{

    DomainSnapshotDataCollection(
        mfem::ParFiniteElementSpace& fes, 
        Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& fields) :
        mesh{*fes.GetMesh()},
        Ex{fields.get(E,X)},
        Ey{fields.get(E,Y)},
        Ez{fields.get(E,Z)},
        Hx{fields.get(H,X)},
        Hy{fields.get(H,Y)},
        Hz{fields.get(H,Z)}
        {}

    void Save(const std::string& folder)
    {
        std::ofstream osex(folder + "/Ex.gf");
        Ex.Save(osex);
        std::ofstream osey(folder + "/Ey.gf");
        Ey.Save(osey);
        std::ofstream osez(folder + "/Ez.gf");
        Ez.Save(osez);
        std::ofstream oshx(folder + "/Hx.gf");
        Hx.Save(oshx);
        std::ofstream oshy(folder + "/Hy.gf");
        Hy.Save(oshy);
        std::ofstream oshz(folder + "/Hz.gf");
        Hz.Save(oshz);
    }

    mfem::Mesh& mesh;
    mfem::ParGridFunction& Ex;
    mfem::ParGridFunction& Ey;
    mfem::ParGridFunction& Ez;
    mfem::ParGridFunction& Hx;
    mfem::ParGridFunction& Hy;
    mfem::ParGridFunction& Hz;
    std::string prefixPath;
};

class FieldProbe {
public:
    FieldProbe(const FieldType& ft, const Direction& d, const Point& p, const int visSteps = 1, const bool writeFile = true) :
        fieldToExtract_{ ft },
        directionToExtract_{ d },
        point_{ p },
        visSteps_{ visSteps },
        write{writeFile}
    {}

    bool write;

    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    const FieldMovie& getFieldMovie() const { return fieldMovie_; }
    const Point& getPoint() const { return point_; }
    int getVisSteps() const { return visSteps_; }
    void addFieldToMovies(double time, const double& field) { fieldMovie_.emplace(time, field); };
    void setProbeID(const size_t id) {id_ = id;}
    size_t getProbeID() const { return id_; }

    std::pair<Time, double> findFrameWithMax() const
    {
        std::pair<Time, double> res{ 0.0, -std::numeric_limits<double>::infinity() };
        for (const auto& [t, f] : fieldMovie_) {
            if (res.second < f) {
                res = { t,f };
            }
        }
        return res;
    }

    std::pair<Time, double> findFrameWithMin() const
    {
        std::pair<Time, double> res{ 0.0, std::numeric_limits<double>::infinity() };
        for (const auto& [t, f] : fieldMovie_) {
            if (res.second > f) {
                res = { t, f };
            }
        }
        return res;
    }
    
private:
    FieldType fieldToExtract_;
    Direction directionToExtract_;
    Point point_;
    int visSteps_;
    size_t id_;

    FieldMovie fieldMovie_;
};

class PointProbe {
public:
    PointProbe(const Point& p, const int visSteps = 1, const bool writeFile = true) :
        point_{ p },
        visSteps_{ visSteps },
        write{writeFile}
    {}

    bool write;

    const FieldMovies& getFieldMovies() const { return fieldMovies_; }
    const Point& getPoint() const { return point_; }
    int getVisSteps() const { return visSteps_; }
    void addFieldsToMovies(Time t, const FieldsForMovie& fields) { fieldMovies_.emplace(t, fields); };
    void setProbeID(const size_t id) {id_ = id;}
    size_t getProbeID() const { return id_; }

private:
    Point point_;
    int visSteps_;
    size_t id_;

    FieldMovies fieldMovies_;
};

struct Probes {
    std::vector<FieldProbe> fieldProbes;
    std::vector<ExporterProbe> exporterProbes;
    std::vector<PointProbe> pointProbes;
    std::vector<NearFieldProbe> nearFieldProbes;
    std::vector<DomainSnapshotProbe> domainSnapshotProbes;
};

}