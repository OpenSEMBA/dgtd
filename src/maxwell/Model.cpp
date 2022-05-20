#include "Model.h"

namespace maxwell {

	Model::Model(Mesh& mesh, const AttributeToMaterial& matMap, const AttributeToBoundary& bdrMap) :
		mesh_(mesh)
	{

		if (matMap.size() == 0) {
			attToMatMap_.emplace(1, Material(1.0, 1.0));
		}
		else {
			attToMatMap_ = matMap;
		}

		if (bdrMap.size() == 0) {
			for (int i = 1; i <= mesh.bdr_attributes.Size(); i++) {
				attToBdrMap_.emplace(i, BdrCond::PEC);
			}
		}
		else {
			attToBdrMap_ = bdrMap;
		}

		if (mesh.bdr_attributes.Size() != attToBdrMap_.size()) {
			throw std::exception("Mesh Boundary Attributes Size and Boundary maps must have same size.");
		}

	}

}
