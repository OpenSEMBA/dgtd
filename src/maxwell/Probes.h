#pragma once

#include <mfem.hpp>
#include "Types.h"

namespace maxwell {

struct FieldViews {
    std::array<mfem::GridFunction*, 3> E;
    std::array<mfem::GridFunction*, 3> H;
};

class Probe {

};

struct ExporterProbe : public Probe {
public:
    ExporterProbe(const std::string& name) :
        name{ name }
    {}

    std::string name{"MaxwellView"};
};

class PointsProbe : public Probe {

public:
    PointsProbe(const FieldType&, const Direction&, const Points&);

    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    const FieldMovie& getFieldMovie() const { return fieldMovie_; }
    const Points& getPoints() const { return points_; }
    void addFrame(double time, const FieldFrame& frame) { fieldMovie_.emplace(time, frame); };
    
private:
    FieldType fieldToExtract_;
    Direction directionToExtract_;
    Points points_;

    FieldMovie fieldMovie_;
};

struct Probes {
    std::vector<PointsProbe> pointsProbes;
    std::vector<ExporterProbe> exporterProbes;

    int visSteps{ 1 };
};

}