#include <components/RCSExporter.h>


namespace maxwell {

using namespace mfem;


RCSCollection::RCSCollection(const NearToFarFieldProbe& p, const SolverOptions& opts) : 
	DataCollection::DataCollection(p.name)
{
}

void RCSCollection::Save()
{
}

}

