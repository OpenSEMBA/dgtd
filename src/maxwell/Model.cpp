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
			for (int i = 0; i < mesh.bdr_attributes.Max(); i++) {
				attToBdrMap_.emplace(i, BdrCond::PEC);
			}
		}
		else {
			attToBdrMap_ = bdrMap;
		}

		if (mesh.bdr_attributes.Max() != attToBdrMap_.size()) {
			throw std::exception("Max Mesh Boundary Attributes and Boundary maps must have same size.");
		}

		bdrMarkers_.SetSize(attToBdrMap_.size());
		bdrMarkers_ = 0; //No condition for any att. Mfem. Obscure.

		for (auto const& imap : attToBdrMap_)
			bdrCondVec_.push_back(imap.second);
	}

}
