#pragma once

#include "components/Sources.h"
#include "evolution/Fields.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, mfem::FiniteElementSpace&);  

    void setInitialFields(Fields&);
    std::array<std::array<mfem::GridFunction, 3>, 2> evalTimeVarField(const Time);
    std::array<std::array<mfem::GridFunction, 3>, 2> evalTimeVarField(const Time, bool is_tf);
    void initTFSFmeshes(const std::pair<SubMesh, SubMesh>& ms);

    Sources sources;

private:

    mfem::FiniteElementSpace& fes_;
    mfem::FiniteElementSpace& tf_fes_,& sf_fes_;
    mfem::SubMesh& tf_mesh_,& sf_mesh_;

};

}