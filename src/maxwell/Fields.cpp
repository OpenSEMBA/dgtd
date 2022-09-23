#include "Fields.h"

#include "MaxwellEvolution3D.h"

namespace maxwell {

using namespace mfem;

Fields::Fields(mfem::FiniteElementSpace& fes) :
    allDOFs{ 
        MaxwellEvolution::numberOfFieldComponents *
        MaxwellEvolution::numberOfMaxDimensions *
        fes.GetNDofs() 
    }
{
    allDOFs = 0.0;

    for (int d = X; d <= Z; d++) {
        E[d].SetSpace(&fes);
        H[d].SetSpace(&fes);
        E[d].SetData(allDOFs.GetData() + d * fes.GetNDofs());
        H[d].SetData(allDOFs.GetData() + (d + 3) * fes.GetNDofs());
    }
}

}