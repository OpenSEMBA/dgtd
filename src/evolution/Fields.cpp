#include "Fields.h"

namespace maxwell {

using namespace mfem;

Fields::Fields(ParFiniteElementSpace& fes)
{
    ParMesh pmesh(MPI_COMM_WORLD, *fes.GetMesh());
    auto fecdg = dynamic_cast<const DG_FECollection*>(fes.FEColl());
    auto fec = DG_FECollection(fes.FEColl()->GetOrder(), pmesh.Dimension(), fecdg->GetBasisType());
    global_fes_ = std::make_unique<ParFiniteElementSpace>(&pmesh, &fec, 3);
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

    e_global_.SetVector(e_[X], 0);
    e_global_.SetVector(e_[Y], fes.GetNDofs());
    e_global_.SetVector(e_[Z], 2 * fes.GetNDofs());
    h_global_.SetVector(h_[X], 0);
    h_global_.SetVector(h_[Y], fes.GetNDofs());
    h_global_.SetVector(h_[Z], 2 * fes.GetNDofs());

}

ParGridFunction& Fields::get(const FieldType& f, const Direction& d)
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

const ParGridFunction& Fields::get(const FieldType& f, const Direction& d) const
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

ParGridFunction& Fields::get(const FieldType& f)
{
    assert(f == E || f == H);
    if (f == E) {
        return e_global_;
    }
    else {
        return h_global_;
    }
}

void Fields::updateGlobal() 
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