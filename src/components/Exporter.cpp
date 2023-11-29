#include "Exporter.h"

namespace maxwell {

using namespace mfem;


NearToFarFieldDataCollection::NearToFarFieldDataCollection(
	const NearToFarFieldProbe& p, const DG_FECollection& fec, FiniteElementSpace& fes, Fields& global)
	: DataCollection(p.name, fes.GetMesh())
{

}


}