#pragma once
#include "mfem.hpp"
#include "Types.h"

namespace maxwell {

class Probe {

public:

    Probe(const FieldType&, const Direction&, DenseMatrix& integPointMat);
    const FieldType& getFieldType() const { return fieldToExtract_; }
    const Direction& getDirection() const { return directionToExtract_; }
    DenseMatrix& getIntegPointMat() { return integPointMat_; }
    FieldMovie& getFieldMovie() { return fieldMovie_; }
    
private:

    FieldType fieldToExtract_;
    Direction directionToExtract_;
    DenseMatrix integPointMat_;
    FieldMovie fieldMovie_;
    
};

struct Probes {
public:

    int vis_steps = 1;
    int precision = 8;
    bool paraview = false;
    bool glvis = false;
    bool extractDataAtPoints = true;
    Probes() = default;

    void addProbeToVector(const Probe& probe) { probeVector_.push_back(probe); }
    std::vector<Probe>& getProbeVector() { return probeVector_; }

private:

    std::vector<Probe> probeVector_;

};

}