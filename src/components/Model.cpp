#include "Model.h"

namespace maxwell {

Model::Model(Mesh& mesh, const GeomTagToMaterial& matMap, const GeomTagToBoundary& bdrMap, const GeomTagToInteriorConditions& intBdrMap) :
	mesh_(mesh)
{
	if (matMap.size() == 0) {
		attToMatMap_.emplace(1, Material(1.0, 1.0));
	}
	else {
		attToMatMap_ = matMap;
	}

	auto f2bdr{ mesh.GetFaceToBdrElMap() };

	for (auto i = bdrMap.begin(); i != bdrMap.end(); i++) {
		faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
	}
	attToBdrMap_ = bdrMap;
	
	if (intBdrMap.size() != 0)
	{
		for (auto i = intBdrMap.begin(); i != intBdrMap.end(); i++){
			if (i->second == BdrCond::PEC || i->second == BdrCond::PMC || i->second == BdrCond::SMA) {
				faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
				attToIntBdrMap_.insert(std::make_pair( i->first, i->second ));
			}
			else {
				throw std::exception("Wrongly declared BdrCond as value in AttributeToInteriorConditions.");
			}
		}
	}

	assembleGeomTagToTypeMap(attToBdrMap_, false);
	assembleGeomTagToTypeMap(attToIntBdrMap_, true);
	assembleBdrToMarkerMaps();

}

void Model::assembleBdrToMarkerMaps()
{
	if (pecMarker_.Size() != 0) {
		bdrToMarkerMap_.insert(std::make_pair(BdrCond::PEC, pecMarker_));
	}
	if (pmcMarker_.Size() != 0) {
		bdrToMarkerMap_.insert(std::make_pair(BdrCond::PMC, pmcMarker_));
	}
	if (smaMarker_.Size() != 0) {
		bdrToMarkerMap_.insert(std::make_pair(BdrCond::SMA, smaMarker_));
	}
	if (intpecMarker_.Size() != 0) {
		intBdrToMarkerMap_.insert(std::make_pair(BdrCond::PEC, intpecMarker_));
	}
	if (intpmcMarker_.Size() != 0) {
		intBdrToMarkerMap_.insert(std::make_pair(BdrCond::PMC, intpmcMarker_));
	}
	if (intsmaMarker_.Size() != 0) {
		intBdrToMarkerMap_.insert(std::make_pair(BdrCond::SMA, intsmaMarker_));
	}
	if (tfsfMarker_.Size() != 0) {
		intBdrToMarkerMap_.insert(std::make_pair(BdrCond::TotalFieldIn, tfsfMarker_));
	}
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

void initMarker(BoundaryMarker& marker, const int size)
{
	marker.SetSize(size);
	marker = 0;
}

void Model::assembleGeomTagToTypeMap(
	std::map<GeomTag, BdrCond>& geomTagToCond,
	bool isInterior)
{
	for (const auto& [geomTag, bdr] : geomTagToCond) {

		if (geomTag <= 0) {
			throw std::exception("geomTag <= 0 in GeomTagToTypeMap assembly.");
		}

		auto& marker{ getMarker(bdr, isInterior) };

		if (marker.Size() == 0) {
			initMarker(getMarker(bdr, isInterior), mesh_.bdr_attributes.Max());
		}

		marker[geomTag - 1] = 1;
	}
}

BoundaryMarker& Model::getMarker(const BdrCond& bdrCond, bool isInterior)
{
	switch (bdrCond) {
	case BdrCond::PEC:
		switch (isInterior) {
			case true:
				return intpecMarker_;
			case false:
				return pecMarker_;
		}
		break;
	case BdrCond::PMC:
		switch (isInterior) {
			case true:
				return intpmcMarker_;
			case false:
				return pmcMarker_;
		}
		break;
	case BdrCond::SMA:
		switch (isInterior) {
			case true:
				return intsmaMarker_;
			case false:
				return smaMarker_;
		}
		break;
	case BdrCond::TotalFieldIn:
		switch (isInterior) {
			case true:
				return tfsfMarker_;
		}
		break;
	default:
		throw std::exception("Wrong BdrCond in getMarkerForBdrCond getter.");
	}
}

}
