#include <components/NearFieldExporter.h>

namespace maxwell {

using namespace mfem;

NearFieldCollection::NearFieldCollection(const SurfaceProbe& p, const SolverOptions& opts) : 
	DataCollection::DataCollection(p.name)
{
}

void NearFieldCollection::Save()
{
}

}

