#pragma once

#include "components/Types.h"

#include <cmath>
#include <memory>

namespace maxwell {

using namespace mfem;

template <typename FES, typename GF>
class Fields {
public:
    Fields(FES& fes, int n_aux = 0);

    GF& get(const FieldType&, const Direction&);
    const GF& get(const FieldType&, const Direction&) const;
    GF& get(const FieldType&);

    mfem::Vector& allDOFs() { return all_dofs_; }
    const mfem::Vector& allDOFs() const { return all_dofs_; }

    int ndofsPerComponent() const { return ndofs_per_component_; }
    int fieldBlockSize() const { return field_block_size_; }
    int auxSize() const { return aux_size_; }

    double getNorml2() const;

private:
    mfem::Vector all_dofs_;
    int ndofs_per_component_ = 0;
    int field_block_size_ = 0;
    int aux_size_ = 0;
    std::unique_ptr<FES> global_fes_;
    std::array<GF, 3> e_, h_;
    GF e_global_, h_global_;
    std::unique_ptr<mfem::DG_FECollection> fec_;

};

template <typename FES, typename GF>
Fields<FES, GF>::Fields(FES& fes, int n_aux)
    : ndofs_per_component_(fes.GetNDofs()),
      field_block_size_(6 * ndofs_per_component_),
      aux_size_(n_aux)
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
    all_dofs_.SetSize(field_block_size_ + aux_size_);
    all_dofs_ = 0.0;
    for (int d = X; d <= Z; d++) {
        e_[d].UseDevice(true);
        h_[d].UseDevice(true);
        e_[d].SetSpace(&fes);
        h_[d].SetSpace(&fes);
        e_[d].MakeRef(all_dofs_,     d  * ndofs_per_component_, ndofs_per_component_);
        h_[d].MakeRef(all_dofs_,(d + 3) * ndofs_per_component_, ndofs_per_component_);
    }

    e_global_.SetSpace(global_fes_.get());
    h_global_.SetSpace(global_fes_.get());

    e_global_.MakeRef(all_dofs_,       0, field_block_size_ / 2);
    h_global_.MakeRef(all_dofs_, field_block_size_ / 2, field_block_size_ / 2);

}

template <typename FES, typename GF>
double Fields<FES, GF>::getNorml2() const
{
    return all_dofs_.Norml2();
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

struct TransferMaps {

    std::array<std::array<mfem::TransferMap, 3>, 2> maps;

    TransferMaps(Fields<ParFiniteElementSpace, ParGridFunction>& src, Fields<FiniteElementSpace, GridFunction>& dst) :
        maps{ std::array<mfem::TransferMap, 3>{ mfem::TransferMap(src.get(E, X), dst.get(E, X)),
                                                 mfem::TransferMap(src.get(E, Y), dst.get(E, Y)),
                                                 mfem::TransferMap(src.get(E, Z), dst.get(E, Z)) },
              std::array<mfem::TransferMap, 3>{ mfem::TransferMap(src.get(H, X), dst.get(H, X)),
                                                 mfem::TransferMap(src.get(H, Y), dst.get(H, Y)),
                                                 mfem::TransferMap(src.get(H, Z), dst.get(H, Z)) } }
    {}

    void transferFields(const Fields<ParFiniteElementSpace, ParGridFunction>& src, Fields<FiniteElementSpace, GridFunction>& dst)
    {
        for (auto f : { E, H }) {
            for (auto d : { X, Y, Z }) {
                maps[f][d].Transfer(src.get(f, d), dst.get(f, d));
            }
        }
    }
};

}