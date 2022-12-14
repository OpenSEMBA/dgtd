#pragma once

#include "Sources.h"
#include "Fields.h"
#include "Types.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(Sources, const mfem::FiniteElementSpace&);  

    void setFields1D(Fields&);
    void setFields3D(Fields&);
    void setGaussianSource(std::unique_ptr<Source> source);

    Sources sources;
private:
    const mfem::FiniteElementSpace& fes_;
};

}