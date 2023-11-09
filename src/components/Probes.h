#pragma once

#include <mfem.hpp>

#include "Types.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "components/SubMesher.h"
#include "components/Model.h"

namespace maxwell {

using FiniteElementOperator = std::unique_ptr<mfemExtension::BilinearForm>;

struct ExporterProbe {
    std::string name{"MaxwellView"};
    int visSteps{ 10 };
};

class PointProbe {
public:
    PointProbe(const FieldType& ft, const Direction& d, const Point& p) :
        fieldToExtract_{ ft },
        directionToExtract_{ d },
        point_{ p }
    {}

    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    const FieldMovie& getFieldMovie() const { return fieldMovie_; }
    const Point& getPoint() const { return point_; }
    void addFieldToMovies(double time, const double& field) { fieldMovie_.emplace(time, field); };

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

    FieldMovie fieldMovie_;
};

class FieldProbe {
public:
    FieldProbe(const Point& p) :
        point_{ p }
    {}

    const PointFieldMovies& getFieldMovies() const { return fieldMovies_; }
    const Point& getPoint() const { return point_; }
    void addFieldsToMovies(Time t, const FieldsForProbes& fields) { fieldMovies_.emplace(t, fields); };


private:

    Point point_;
    PointFieldMovies fieldMovies_;

};

class NearToFarFieldProbe {
public:

    NearToFarFieldProbe();
    NearToFarFieldProbe(const mfem::Array<int>& tags);

    const Array<int>& getTags() { return tags_; }

    void buildSubMesher(const Mesh& mesh, const Array<int>& marker);
    NearToFarFieldSubMesher* getSubMesher() { return ntff_sm_.get(); }

    void addFieldsToMovies(Time t, const GridFuncForProbes& fields) { faceFieldMovies_.emplace(t, fields); };
    const FaceFieldMovies& getFieldMovies() { return faceFieldMovies_; }

private:

    FaceFieldMovies faceFieldMovies_;
    Array<int> tags_;
    std::unique_ptr<NearToFarFieldSubMesher> ntff_sm_;

};

struct Probes {
    std::vector<PointProbe>          pointProbes;
    std::vector<ExporterProbe>       exporterProbes;
    std::vector<FieldProbe>          fieldProbes;
    std::vector<NearToFarFieldProbe> nearToFarFieldProbes;
};

}