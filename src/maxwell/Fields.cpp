#include "Fields.h"
#include "evolution/Evolution3D.h"

namespace maxwell {

using namespace mfem;

Fields::Fields(mfem::FiniteElementSpace& fes)
{
    allDOFs_.SetSize(MaxwellEvolution3D::numberOfFieldComponents * MaxwellEvolution3D::numberOfMaxDimensions * fes.GetNDofs());
    allDOFs_ = 0.0;
    for (int d = X; d <= Z; d++) {
        e_[d].SetSpace(&fes);
        h_[d].SetSpace(&fes);
        e_[d].SetDataAndSize(allDOFs_.GetData() + d * fes.GetNDofs(),       fes.GetNDofs());
        h_[d].SetDataAndSize(allDOFs_.GetData() + (d + 3) * fes.GetNDofs(), fes.GetNDofs());
    }
}

GridFunction& Fields::operator()(const FieldType& f, const Direction& d)
{
    return get(f, d);
}

GridFunction& Fields::get(const FieldType& f, const Direction& d)
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

}