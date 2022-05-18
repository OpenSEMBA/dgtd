#include "Model.h"

namespace maxwell {

Model::Model(Mesh& mesh, const AttributeToMaterial& matVec) :
mesh_(mesh),
attToMatVec_(matVec)
{
}

}
