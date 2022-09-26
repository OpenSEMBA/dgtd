#include "Fields.h"
#include "MaxwellEvolution1D.h"
#include "MaxwellEvolution3D.h"

namespace maxwell {

using namespace mfem;

Fields::Fields(mfem::FiniteElementSpace& fes)
{
    allDOFs = 0.0;
    
    switch(fes.GetMesh()->Dimension()) {
    case 1:
        allDOFs = MaxwellEvolution1D::numberOfFieldComponents * MaxwellEvolution1D::numberOfMaxDimensions * fes.GetNDofs();
        break;
    default:
        allDOFs = MaxwellEvolution3D::numberOfFieldComponents * MaxwellEvolution3D::numberOfMaxDimensions * fes.GetNDofs();
    }

    for (int d = X; d <= Z; d++) {
        E[d].SetSpace(&fes);
        H[d].SetSpace(&fes);
        E[d].SetData(allDOFs.GetData() + d * fes.GetNDofs());
        H[d].SetData(allDOFs.GetData() + (d + 3) * fes.GetNDofs());
    }
}

}