#include "Model.h"

namespace maxwell {

	Model::Model(Mesh& mesh, const AttributeToMaterial& matMap, const AttributeToBoundary& bdrMap) :
		mesh_(mesh)
	{
		if (matMap.size() != bdrMap.size()) {
			throw std::exception("Material and Boundary maps must have same size.");
		}

		if (matMap.size() == 0) {
			attToMatMap_.emplace(1, Material(1.0, 1.0));
		}
		else {
			attToMatMap_ = matMap;
		}

		if (bdrMap.size() == 0) {
			attToBdrMap_.emplace(1, BdrCond::PEC);
		}
		else {
			attToBdrMap_ = bdrMap;
		}
		for (auto const& it : attToBdrMap_) {
			bdrMarkers_.Append(it.first);
			bdrCondArr_.Append(it.second);
		}
}

}
