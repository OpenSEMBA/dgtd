#include "Probes.h"

namespace maxwell {

void NearToFarFieldProbe::buildSubMesher(const Mesh& mesh, const Array<int>& marker)
{
	auto sm{ NearToFarFieldSubMesher(mesh, marker) };
	ntff_sm_ = std::make_unique<NearToFarFieldSubMesher>(sm);
}

}