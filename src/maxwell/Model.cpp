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

	for (const auto& kv : attToBdrMap_) {
		const auto& att{ kv.first };
		const auto& bdr{ kv.second };
		assert(att > 0);

		BoundaryMarker bdrMarker{ mesh_.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[att - 1] = 1;
		
		bdrToMarkerMap_.emplace(bdr, bdrMarker);
	}
}


mfem::Vector Model::buildPiecewiseArgVector(const FieldType& f) const
{
	mfem::Vector res;
	res.SetSize((int)attToMatMap_.size());

	int i = 0;
	for (auto const& kv : attToMatMap_) {
		switch (f) {
		case FieldType::E:
			res[i] = kv.second.getPermittivity();
			break;
		case FieldType::H:
			res[i] = kv.second.getPermeability();
			break;
		}
		i++;
	}

	return res;
}


}
