#pragma once

#include <mfem.hpp>
#include "Types.h"

namespace maxwell {

struct FieldViews {
    std::array<mfem::GridFunction, 3>* E;
    std::array<mfem::GridFunction, 3>* H;
};

class Probe {

};

class ExporterProbe : public Probe {
};

class PointsProbe : public Probe {

public:
    PointsProbe(const FieldType&, const Direction&, std::vector<std::vector<double>>& integPoints);
    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    const mfem::DenseMatrix& getIntegPointMat() const { return integPointMat_; }
    const FieldMovie& getFieldMovie() const { return fieldMovie_; }

    void addFrame(double time, const FieldFrame& frame) { fieldMovie_.emplace(time, frame); };
    
private:
    FieldType fieldToExtract_;
    Direction directionToExtract_;
    mfem::DenseMatrix integPointMat_;
    FieldMovie fieldMovie_;

    const bool verifyEntryVectorsSameSize(std::vector<std::vector<double>>& points) const;
    const void verifyEntrySubvectorsNotEmpty(std::vector<std::vector<double>>& points) const;
    const void buildIntegPointMat(std::vector<std::vector<double>>& points);
};

struct Probes {
    std::vector<PointsProbe> pointsProbes;
    std::vector<ExporterProbe> exporterProbes;
};

}