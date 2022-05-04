#include "Model.h"
#include <stdexcept>

namespace maxwell {

Model::Model(Mesh mesh, std::map<attribute, Material> matMap) :
mesh_(mesh),
attToMatMap_(matMap)
{
}

}
