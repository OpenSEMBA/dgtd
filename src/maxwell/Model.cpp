#include "Model.h"
#include <stdexcept>

namespace maxwell {

Model::Model(Mesh mesh, attToMaterialMap matMap) :
mesh_(mesh),
matMap_(matMap)
{
	if (mesh_.Dimension() != 1)
		throw std::runtime_error("Mesh dimension is not one.");

}

}
