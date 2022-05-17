#pragma once
#include "mfem.hpp"
#include "Types.h"

namespace maxwell {

class Probe {
public:

    Probe(const FieldType&, const Direction&, const DenseMatrix&);
    FieldType& getFieldType() { return fieldToExtract_; }
    Direction& getDirection() { return directionToExtract_; }
    DenseMatrix& getIntegPointMat() { return integPointMat_; }
private:
    FieldType fieldToExtract_;
    Direction directionToExtract_;
    DenseMatrix integPointMat_;
    
};

struct Probes {
public:

    int vis_steps = 1;
    int precision = 8;
    bool paraview = false;
    bool glvis = false;
    bool extractDataAtPoints_ = true;
    Probes() = default;

    void addProbeToVector(const Probes& source) { probesVector_.push_back(source); }
    std::vector<Probes> getProbesVector() const { return probesVector_; }

private:

    std::vector<Probes> probesVector_;
};

}