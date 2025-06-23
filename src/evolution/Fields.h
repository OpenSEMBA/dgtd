#pragma once

#include <mfem.hpp>

#include "components/Types.h"

namespace maxwell {

using namespace mfem;

template <typename FES, typename GF>
class Fields {
public:
    Fields(FES& fes);
    
    GF& get(const FieldType&, const Direction&);
    const GF& get(const FieldType&, const Direction&) const;
    GF& get(const FieldType&);
    
    mfem::Vector& allDOFs() { return allDOFs_; }
    const mfem::Vector& allDOFs() const { return allDOFs_; }

    double getNorml2() const { return allDOFs_.Norml2(); }

    void updateGlobal();

private:
    mfem::Vector allDOFs_;
    std::unique_ptr<FES> global_fes_;
    std::array<GF, 3> e_, h_;
    GF e_global_, h_global_;
    std::unique_ptr<mfem::DG_FECollection> fec_;

};

template <typename FES, typename GF>
Fields<FES, GF>::Fields(FES& fes)
{
    auto fecdg = dynamic_cast<const DG_FECollection*>(fes.FEColl());
    fec_ = std::make_unique<DG_FECollection>(fes.FEColl()->GetOrder(), fes.GetMesh()->Dimension(), fecdg->GetBasisType());
        if constexpr (std::is_same<FES, FiniteElementSpace>::value)
    {
        global_fes_ = std::make_unique<FiniteElementSpace>(fes.GetMesh(), fec_.get(), 3);
    }
    else if constexpr (std::is_same<FES, ParFiniteElementSpace>::value)
    {
         global_fes_ = std::make_unique<ParFiniteElementSpace>(fes.GetParMesh(), fec_.get(), 3);
    }
    allDOFs_.SetSize(6 * fes.GetNDofs());
    allDOFs_ = 0.0;
    for (int d = X; d <= Z; d++) {
        e_[d].SetDataAndSize(allDOFs_.GetData() + d *       fes.GetNDofs(), fes.GetNDofs());
        h_[d].SetDataAndSize(allDOFs_.GetData() + (d + 3) * fes.GetNDofs(), fes.GetNDofs());
        e_[d].SetSpace(&fes);
        h_[d].SetSpace(&fes);
    }

    e_global_.SetSpace(global_fes_.get());
    h_global_.SetSpace(global_fes_.get());

    auto dofsize = global_fes_->GetNDofs() / 3;

    e_global_.SetVector(e_[X], 0);
    e_global_.SetVector(e_[Y], dofsize);
    e_global_.SetVector(e_[Z], 2 * dofsize);
    h_global_.SetVector(h_[X], 0);
    h_global_.SetVector(h_[Y], dofsize);
    h_global_.SetVector(h_[Z], 2 * dofsize);

    
}

template <typename FES, typename GF>
GF& Fields<FES, GF>::get(const FieldType& f, const Direction& d)
{
    assert(f == E || f == H);
    assert(d == X || d == Y || d == Z);
    if (f == E) {
        return e_[d];
    }
    else {
        return h_[d];
    }
}

template <typename FES, typename GF>
const GF& Fields<FES, GF>::get(const FieldType& f, const Direction& d) const
{
    assert(f == E || f == H);
    assert(d == X || d == Y || d == Z);
    if (f == E) {
        return e_[d];
    }
    else {
        return h_[d];
    }
}

template <typename FES, typename GF>
GF& Fields<FES, GF>::get(const FieldType& f)
{
    assert(f == E || f == H);
    if (f == E) {
        return e_global_;
    }
    else {
        return h_global_;
    }
}

template <typename FES, typename GF>
void Fields<FES, GF>::updateGlobal() 
{
    auto offset = global_fes_->GetNDofs();
    e_global_.SetVector(e_[X], 0);
    e_global_.SetVector(e_[Y], offset);
    e_global_.SetVector(e_[Z], 2 * offset);
    h_global_.SetVector(h_[X], 0);
    h_global_.SetVector(h_[Y], offset);
    h_global_.SetVector(h_[Z], 2 * offset);
}


}