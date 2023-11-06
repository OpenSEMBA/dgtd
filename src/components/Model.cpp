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
				std::runtime_error("Wrongly declared BdrCond as value in AttributeToInteriorConditions.");
			}
		}
	}

	initGeomTagToTypeMaps();
	assembleGeomTagToTypeMap(attToBdrMap_, false);
	assembleGeomTagToTypeMap(attToIntBdrMap_, true);

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


void Model::assembleGeomTagToTypeMap(
	std::map<GeomTag, BdrCond>& geomTagToCond,
	bool isInterior)
{
	for (const auto& [geomTag, bdr] : geomTagToCond) {

		if (geomTag <= 0) {
			throw std::exception("geomTag <= 0 in GeomTagToTypeMap assembly.");
		}

		switch (isInterior) {
		case false:
			getMarkerForBdrCond(bdr)[geomTag - 1] = 1;
			break;
		case true:
			getInteriorMarkerForBdrCond(bdr)[geomTag - 1] = 1;
			break;
		}
	}
}

void Model::initGeomTagToTypeMaps()
{
	pecMarker_.SetSize(mesh_.bdr_attributes.Max());
	pmcMarker_.SetSize(mesh_.bdr_attributes.Max());
	smaMarker_.SetSize(mesh_.bdr_attributes.Max());
	intpecMarker_.SetSize(mesh_.bdr_attributes.Max());
	intpmcMarker_.SetSize(mesh_.bdr_attributes.Max());
	intsmaMarker_.SetSize(mesh_.bdr_attributes.Max());
	
	pecMarker_    = 0;
	pmcMarker_    = 0;
	smaMarker_    = 0;
	intpecMarker_ = 0;
	intpmcMarker_ = 0;
	intsmaMarker_ = 0;

}

BoundaryMarker& Model::getMarkerForBdrCond(const BdrCond& bdrCond)
{
	switch (bdrCond) {
	case BdrCond::PEC:
		return pecMarker_;
	case BdrCond::PMC:
		return pmcMarker_;
	case BdrCond::SMA:
		return smaMarker_;
	default:
		throw std::exception("Wrong BdrCond in getMarkerForBdrCond getter.");
	}
}

InteriorBoundaryMarker& Model::getInteriorMarkerForBdrCond(const BdrCond& bdrCond)
{
	switch (bdrCond) {
	case BdrCond::PEC:
		return intpecMarker_;
	case BdrCond::PMC:
		return intpmcMarker_;
	case BdrCond::SMA:
		return intsmaMarker_;
	default:
		throw std::exception("Wrong BdrCond in getInteriorMarkerForBdrCond getter.");
	}
}


}
