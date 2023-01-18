#pragma once

#include <mfem.hpp>

#include "Types.h"

namespace maxwell {

struct ExporterProbe {
public:
    ExporterProbe(const std::string& name) :
        name{ name }
    {}

    std::string name{"MaxwellView"};
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
    void addFieldToMovie(double time, const double& field) { fieldMovie_.emplace(time, field); };
    
private:
    FieldType fieldToExtract_;
    Direction directionToExtract_;
    Point point_;

    FieldMovie fieldMovie_;
};

struct Probes {
    std::vector<PointProbe> pointProbes;
    std::vector<ExporterProbe> exporterProbes;

    int visSteps{ 10 };
};

}