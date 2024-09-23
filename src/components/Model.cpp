#include "Model.h"

namespace maxwell {

Model::Model(Mesh& mesh, const GeomTagToMaterialInfo& matInfo, const GeomTagToBoundaryInfo& bdrInfo) :
	mesh_(mesh)
{
	if (matInfo.gt2m.size() == 0) {
		attToMatMap_.emplace(1, Material(1.0, 1.0, 0.0));
	}
	else {
		attToMatMap_ = matInfo.gt2m;
	}

	auto f2bdr{ mesh.GetFaceToBdrElMap() };

	for (auto i = bdrInfo.gt2b.begin(); i != bdrInfo.gt2b.end(); i++) {
		faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
	}
	attToBdrMap_ = bdrInfo.gt2b;
	
	if (bdrInfo.gt2ib.size() != 0)
	{
		for (auto i = bdrInfo.gt2ib.begin(); i != bdrInfo.gt2ib.end(); i++){
			if (i->second == BdrCond::PEC || i->second == BdrCond::PMC || i->second == BdrCond::SMA) {
				faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
				attToIntBdrMap_.insert(std::make_pair( i->first, i->second ));
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
}

std::size_t Model::numberOfMaterials() const
{
	return attToMatMap_.size();
}

std::size_t Model::numberOfBoundaryMaterials() const
{
	return attToBdrMap_.size();
}

mfem::Vector Model::initialiseGeomTagVector() const
{
	int size{ 0 };
	for (auto const& [geomTag, mat] : attToMatMap_) {
		if (geomTag >= size) {
			size = geomTag;
		}
	}
	assert(size > 0);
	mfem::Vector res(size);
	res = 0.0;
	return res;
}

mfem::Vector Model::buildEpsMuPiecewiseVector(const FieldType& f) const
{
	auto res{ initialiseGeomTagVector() };

	for (auto const& [geomTag, mat] : attToMatMap_) {
		switch (f) {
		case FieldType::E:
			res[geomTag - 1] = mat.getPermittivity();
			break;
		case FieldType::H:
			res[geomTag - 1] = mat.getPermeability();
			break;
		}
	}

	return res;
}

mfem::Vector Model::buildSigmaPiecewiseVector() const
{
	auto res{ initialiseGeomTagVector() };

	for (auto const& [geomTag, mat] : attToMatMap_) {
		res[geomTag - 1] = mat.getConductivity();
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
			throw std::runtime_error("geomTag <= 0 in GeomTagToTypeMap assembly.");
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
			return tfsfMarker_;
		break;
	default:
		throw std::runtime_error("Wrong BdrCond in getMarkerForBdrCond getter.");
	}
}

}
