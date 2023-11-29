#pragma once

#include <evolution/Fields.h>
#include <components/Types.h>
#include <components/SubMesher.h>
#include <components/Probes.h>

namespace maxwell {

using namespace mfem;

class NearToFarFieldDataCollection : public DataCollection
{
public:

	NearToFarFieldDataCollection(const NearToFarFieldProbe&, const DG_FECollection& fec, FiniteElementSpace& fes, Fields&);
	
};

}