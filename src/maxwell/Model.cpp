#include "Model.h"

namespace maxwell {

Model::Model(Mesh mesh, std::map<attribute, Material> matMap) :
mesh_(mesh),
attToMatMap_(matMap)
{
}

}
