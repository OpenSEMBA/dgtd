#include "Fields.h"
#include "evolution/Evolution3D.h"

namespace maxwell {

using namespace mfem;

Fields::Fields(mfem::FiniteElementSpace& fes)
{
    allDOFs.SetSize(MaxwellEvolution3D::numberOfFieldComponents * MaxwellEvolution3D::numberOfMaxDimensions * fes.GetNDofs());
    allDOFs = 0.0;
    for (int d = X; d <= Z; d++) {
        E[d].SetSpace(&fes);
        H[d].SetSpace(&fes);
        E[d].SetDataAndSize(allDOFs.GetData() + d * fes.GetNDofs(),       fes.GetNDofs());
        H[d].SetDataAndSize(allDOFs.GetData() + (d + 3) * fes.GetNDofs(), fes.GetNDofs());
    }
}

}