#pragma once

#include "Sources.h"
#include "Fields.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(Sources, const mfem::FiniteElementSpace&);  

    void setFields(Fields&);

    Sources sources;
private:
    const mfem::FiniteElementSpace& fes_;
};

}