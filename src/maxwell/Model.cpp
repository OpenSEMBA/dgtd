#include "Model.h"

namespace maxwell {

Model::Model(Mesh& mesh, std::vector<std::pair<attribute, Material>>& matVec) :
mesh_(mesh),
attToMatVec_(matVec)
{
}

}
