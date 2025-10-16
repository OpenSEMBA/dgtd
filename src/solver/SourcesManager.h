#pragma once

#include "components/Sources.h"
#include "evolution/Fields.h"
#include "components/SubMesher.h"

using FieldGridFuncs = std::array<std::array<mfem::GridFunction, 3>, 2>;

namespace maxwell {

class SourcesManager {
public:
    SourcesManager(const Sources&, mfem::ParFiniteElementSpace&, Fields<ParFiniteElementSpace, ParGridFunction>& fields);

    FieldGridFuncs evalTimeVarField(const Time, FiniteElementSpace*);
    void initTFSFPreReqs(const ParMesh&, const Array<int>& marker);
    FiniteElementSpace* getTFSpace() { return tf_fes_.get(); }
    FiniteElementSpace* getSFSpace() { return sf_fes_.get(); }
    FiniteElementSpace* getGlobalTFSFSpace() { return global_tfsf_fes_.get(); }
    TotalFieldScatteredFieldSubMesher& getTFSFSubMesher() { return tfsf_submesher_; }
    void markDoFSforTForSF(FieldGridFuncs&, bool isTF);

    Sources sources;

private:

    void initTFSFSubMesher(const ParMesh&, const Array<int>& marker);
    void initTFSFSpaces();
    void setInitialFields(Fields<ParFiniteElementSpace, ParGridFunction>&);

    mfem::ParFiniteElementSpace& fes_;
    TotalFieldScatteredFieldSubMesher tfsf_submesher_;
    std::unique_ptr<FiniteElementSpace> tf_fes_, sf_fes_, global_tfsf_fes_;
    

};

}