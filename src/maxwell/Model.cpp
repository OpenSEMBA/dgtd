#include "Model.h"

namespace maxwell {

Model::Model(Mesh& mesh, const AttributeToMaterial& matVec, const AttributeToBoundary& bdrVec) :
mesh_(mesh)
{
	if (matVec.size() == 0) {
		attToMatVec_.emplace(1, Material(1.0, 1.0));
	}
	else {
		attToMatVec_ = matVec;
	}

	if (bdrVec.size() == 0) {
		attToBdrVec_.emplace(1, BdrCond::PEC);
	}
	else {
		attToBdrVec_ = bdrVec;
	}
}

}
