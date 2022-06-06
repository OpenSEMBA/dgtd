#pragma once
#include "mfem.hpp"
#include "Types.h"


namespace maxwell {

class Probe {

};

class ExporterProbe : public Probe {
public:
    enum class Type {
        Paraview,
        Glvis
    };

    Type type = Type::Paraview;
    int precision = 8;
};

class PointsProbe : public Probe {

public:

    Probe(const FieldType&, const Direction&, std::vector<std::vector<double>>& integPoints);
    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    DenseMatrix& getIntegPointMat() { return integPointMat_; }
    FieldMovie& getFieldMovie() { return fieldMovie_; }
    const FieldMovie& getConstFieldMovie() const { return fieldMovie_; }
    
private:

    FieldType fieldToExtract_;
    Direction directionToExtract_;
    DenseMatrix integPointMat_;
    FieldMovie fieldMovie_;

    const bool verifyEntryVectorsSameSize(std::vector<std::vector<double>>& points) const;
    const void verifyEntrySubvectorsNotEmpty(std::vector<std::vector<double>>& points) const;
    const void buildIntegPointMat(std::vector<std::vector<double>>& points);
    
};

class Probes {
public:

    int vis_steps = 1;
    bool extractDataAtPoints = false;

    void addProbeToVector(const Probe& probe) { probeVector_.push_back(probe); }
    
    std::vector<PointsProbe&> getPointsProbes();
    std::vector<ExporterProbe&> getExporterProbes();

private:
    boost::base_collection<Probe> probes_;
};

}