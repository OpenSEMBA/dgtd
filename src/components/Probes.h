#pragma once

#include <mfem.hpp>

#include "Types.h"

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
    std::string exportPath = std::string("./defaultExportPath/" + name + "/");
    int expSteps { 10 };
};

class FieldProbe {
public:
    FieldProbe(const FieldType& ft, const Direction& d, const Point& p, const bool writeFile = false) :
        fieldToExtract_{ ft },
        directionToExtract_{ d },
        point_{ p },
        write{writeFile}
    {}

    bool write;

    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    const FieldMovie& getFieldMovie() const { return fieldMovie_; }
    const Point& getPoint() const { return point_; }
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
    size_t id_;

    FieldMovie fieldMovie_;
};

class PointProbe {
public:
    PointProbe(const Point& p, const bool writeFile = false) :
        point_{ p },
        write{writeFile}
    {}

    bool write;

    const FieldMovies& getFieldMovies() const { return fieldMovies_; }
    const Point& getPoint() const { return point_; }
    void addFieldsToMovies(Time t, const FieldsForMovie& fields) { fieldMovies_.emplace(t, fields); };
    void setProbeID(const size_t id) {id_ = id;}
    size_t getProbeID() const { return id_; }

private:

    Point point_;
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