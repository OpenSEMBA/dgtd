#pragma once

#include "components/Sources.h"
#include "evolution/Fields.h"
#include "components/SubMesher.h"

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, mfem::FiniteElementSpace&);  

    void setInitialFields(Fields&);
    std::array<std::array<mfem::GridFunction, 3>, 2> evalTimeVarField(const Time);
    std::array<std::array<mfem::GridFunction, 3>, 2> evalTimeVarField(const Time, bool is_tf);
    void initTFSFPreReqs(const Mesh&);
    FiniteElementSpace* getTFSpace() { return tf_fes_.get(); }
    FiniteElementSpace* getSFSpace() { return sf_fes_.get(); }
    TotalFieldScatteredFieldSubMesher& getTFSFSubMesher() { return tfsf_submesher_; }

    Sources sources;

private:

    void initTFSFSubMesher(const Mesh&);
    void initTFSFSpaces();

    mfem::FiniteElementSpace& fes_;
    TotalFieldScatteredFieldSubMesher tfsf_submesher_;
    std::unique_ptr<FiniteElementSpace> tf_fes_, sf_fes_;
    

};

}