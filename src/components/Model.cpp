#include "Model.h"

namespace maxwell {

Model::Model(Mesh& mesh, const AttributeToMaterial& matMap, const AttributeToBoundary& bdrMap, const AttributeToInteriorConditions& intBdrMap) :
	mesh_(mesh)
{
	if (matMap.size() == 0) {
		attToMatMap_.emplace(1, Material(1.0, 1.0));
	}
	else {
		attToMatMap_ = matMap;
	}

	attToBdrMap_ = bdrMap;
	
	if (intBdrMap.size() != 0)
	{
		for (auto i = intBdrMap.begin(); i != intBdrMap.end(); i++){
			if (i->second == BdrCond::PEC || i->second == BdrCond::PMC || i->second == BdrCond::SMA || i->second == BdrCond::NONE) {
				attToIntBdrMap_.insert({ i->first, i->second });
			}
			else {
				std::runtime_error("Wrongly declared BdrCond as value in AttributeToInteriorConditions");
			}
		}		
	}

	assembleAttToTypeMap(attToBdrMap_, bdrToMarkerMap_);
	assembleAttToTypeMap(attToIntBdrMap_, intBdrToMarkerMap_);
	assembleAttToTypeMap(attToIntSrcMap_, intSrcToMarkerMap_);
}

std::size_t Model::numberOfMaterials() const
{
	return attToMatMap_.size();
}

std::size_t Model::numberOfBoundaryMaterials() const
{
	return attToBdrMap_.size();
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


void Model::assembleAttToTypeMap(
	std::map<Attribute, BdrCond>& attToCond,
	std::multimap<BdrCond, BoundaryMarker>& attToMarker)
{
	for (const auto& kv : attToCond) {
		const auto& att{ kv.first };
		const auto& bdr{ kv.second };
		assert(att > 0);

		BoundaryMarker bdrMarker{ mesh_.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[att - 1] = 1;

		attToMarker.emplace(bdr, bdrMarker);
	}
}


}
