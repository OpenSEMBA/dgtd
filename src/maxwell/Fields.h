#pragma once

#include <mfem.hpp>

#include "Types.h"
#include "Sources.h"

namespace maxwell {

class Fields {
public:
    Fields(mfem::FiniteElementSpace& fes);
    
    std::array<mfem::GridFunction, 3> E, H;
    mfem::GridFunction E1D, H1D;
    mfem::Vector allDOFs;

    double getNorml2() const { return allDOFs.Norml2(); }

};
}