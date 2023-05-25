#pragma once

#include "components/Sources.h"
#include "evolution/Fields.h"

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