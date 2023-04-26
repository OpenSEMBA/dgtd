#pragma once

#include "Sources.h"
#include "Fields.h"
#include "Types.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, mfem::FiniteElementSpace&);  

    void setInitialFields(Fields&);
    std::array<std::array<mfem::GridFunction, 3>, 2> evalTimeVarField(const double time);

    Sources sources;

private:
    mfem::FiniteElementSpace& fes_;
};

}