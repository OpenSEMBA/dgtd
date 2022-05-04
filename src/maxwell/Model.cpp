#include "Model.h"
#include <stdexcept>

namespace maxwell {

Model::Model(Mesh mesh, std::map<attribute, Material> matMap) :
mesh_(mesh),
attToMatMap_(matMap)
{
	if (mesh_.Dimension() != 1)
		throw std::runtime_error("Mesh dimension is not one.");

}

}
