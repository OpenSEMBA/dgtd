#pragma once

#include "Sources.h"
#include "Fields.h"
#include "Types.h"

namespace maxwell {

class SourcesManager {
public:

    using TimeVarOperators = std::array<mfem::GridFunction, 3>;

    SourcesManager(const Sources&, mfem::FiniteElementSpace&);  

    void setInitialFields(Fields&);
    TimeVarOperators evalTimeVarField(const double time);

    Sources sources;

private:
    mfem::FiniteElementSpace& fes_;
};

}