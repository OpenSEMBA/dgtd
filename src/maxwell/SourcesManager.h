#pragma once

#include "Sources.h"
#include "Fields.h"
#include "Types.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, mfem::FiniteElementSpace&);  

    void setInitialFields(Fields&);
    mfem::GridFunction evalTotalField(const double time);

    Sources sources;

private:
    mfem::FiniteElementSpace& fes_;
};

}