#pragma once

#include "Sources.h"
#include "Fields.h"
#include "Types.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, const mfem::FiniteElementSpace&);  

    void setFields1D(Fields&);
    void setFields3D(Fields&);
    
    Sources sources;
private:
    const mfem::FiniteElementSpace& fes_;
};

}