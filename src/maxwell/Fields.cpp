#include "Fields.h"
#include "Evolution1D.h"
#include "Evolution3D.h"

namespace maxwell {

using namespace mfem;

Fields::Fields(mfem::FiniteElementSpace& fes)
{
    switch (fes.GetMesh()->Dimension()) {
    case 1:
        allDOFs.SetSize(MaxwellEvolution1D::numberOfFieldComponents * MaxwellEvolution1D::numberOfMaxDimensions * fes.GetNDofs());
        allDOFs = 0.0;
        E1D.SetSpace(&fes);
        H1D.SetSpace(&fes);
        E1D.SetDataAndSize(allDOFs.GetData(),                  fes.GetNDofs());
        H1D.SetDataAndSize(allDOFs.GetData() + fes.GetNDofs(), fes.GetNDofs());
        break;
    default:
        allDOFs.SetSize(MaxwellEvolution3D::numberOfFieldComponents * MaxwellEvolution3D::numberOfMaxDimensions * fes.GetNDofs());
        allDOFs = 0.0;
        for (int d = X; d <= Z; d++) {
            E[d].SetSpace(&fes);
            H[d].SetSpace(&fes);
            E[d].SetDataAndSize(allDOFs.GetData() + d * fes.GetNDofs(),       fes.GetNDofs());
            H[d].SetDataAndSize(allDOFs.GetData() + (d + 3) * fes.GetNDofs(), fes.GetNDofs());
        }
        break;
    }
}

}