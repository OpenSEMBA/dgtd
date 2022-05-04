#pragma once
#include "mfem.hpp"
#include "Types.h"

namespace maxwell {

class Probes {
    public:
        int vis_steps = 1;
        int precision = 8;
        bool paraview = false;
        bool glvis = false;
        bool extractDataAtPoints = false;
        FieldType fieldToExtract = FieldType::E;
        DenseMatrix integPointMat;
    private:
};

}