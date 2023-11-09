#pragma once

#include <mfem.hpp>

#include "Types.h"
#include "mfemExtension/BilinearIntegrators.h"

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

    const FieldMovies& getFieldMovies() const { return fieldMovies_; }
    const Point& getPoint() const { return point_; }
    void addFieldsToMovies(Time t, const FieldsForFP& fields) { fieldMovies_.emplace(t, fields); };


private:

    Point point_;

    FieldMovies fieldMovies_;
};

//class EnergyProbe {
//public:
//
//    EnergyProbe();
//
//    void addFieldsToMovie(Time t, const Fields& fields) { fieldsMovie_.emplace(t, fields); };
//    double getEnergy(const FieldType& ft, const Direction& d);
//
//    int visSteps{ 10 };
//
//private:
//
//    std::map<Time, Fields> fieldsMovie_;

//};

struct Probes {
    std::vector<PointProbe> pointProbes;
    std::vector<ExporterProbe> exporterProbes;
    std::vector<FieldProbe> fieldProbes;
    //std::vector<EnergyProbe> energyProbes;
};

}