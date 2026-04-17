#pragma once

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
    
    mfem::Vector& allDOFs() { return all_dofs_; }
    const mfem::Vector& allDOFs() const { return all_dofs_; }

    double getNorml2() const { return all_dofs_.Norml2(); }

private:
    mfem::Vector all_dofs_;
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
    all_dofs_.UseDevice(true);
    e_global_.UseDevice(true);
    h_global_.UseDevice(true);
    all_dofs_.SetSize(6 * fes.GetNDofs());
    all_dofs_ = 0.0;
    for (int d = X; d <= Z; d++) {
        e_[d].UseDevice(true);
        h_[d].UseDevice(true);
        e_[d].SetSpace(&fes);
        h_[d].SetSpace(&fes);
        e_[d].MakeRef(all_dofs_,     d  * fes.GetNDofs(), fes.GetNDofs());
        h_[d].MakeRef(all_dofs_,(d + 3) * fes.GetNDofs(), fes.GetNDofs());
    }

    e_global_.SetSpace(global_fes_.get());
    h_global_.SetSpace(global_fes_.get());

    auto field_dof_size = all_dofs_.Size() / 2;

    e_global_.MakeRef(all_dofs_,       0, field_dof_size);
    h_global_.MakeRef(all_dofs_, field_dof_size, field_dof_size);

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

}