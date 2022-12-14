#include "Fields.h"
#include "MaxwellEvolution1D.h"
#include "MaxwellEvolution3D.h"

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
        E1D.SetData(allDOFs.GetData());
        H1D.SetData(allDOFs.GetData() + fes.GetNDofs());
        break;
    default:
        allDOFs.SetSize(MaxwellEvolution3D::numberOfFieldComponents * MaxwellEvolution3D::numberOfMaxDimensions * fes.GetNDofs());
        allDOFs = 0.0;
        for (int d = X; d <= Z; d++) {
            E[d].SetSpace(&fes);
            H[d].SetSpace(&fes);
            E[d].SetData(allDOFs.GetData() + d * fes.GetNDofs());
            H[d].SetData(allDOFs.GetData() + (d + 3) * fes.GetNDofs());
        }
        break;
    }
}

}