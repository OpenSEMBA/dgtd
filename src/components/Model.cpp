#include "Model.h"

namespace maxwell {

using namespace mfem;


std::map<GlobalElementId, Position> buildSerialElem2CenterMap(Mesh& mesh){
	std::map<GlobalElementId, Position> res;
	for (auto e = 0; e < mesh.GetNE(); e++){
		Vector center;
		mesh.GetElementCenter(e, center);
		res[e] = center;
	}
	return res;
}


std::map<LocalElementId, Position> buildPartitionElem2CenterMap(ParMesh& pmesh){
	std::map<LocalElementId, Position> res;
	for (auto e = 0; e < pmesh.GetNE(); e++){
		Vector center;
		pmesh.GetElementCenter(e, center);
		res[e] = center;
	}
	return res;
}

std::map<GlobalElementId, LocalElementId> buildGlobalToPartitionLocalElementMap(const std::map<GlobalElementId, Position>& serial, const std::map<LocalElementId, Position>& local)
{
	double tol = 1e-5;
	std::map<GlobalElementId, LocalElementId> res;
	for (const auto& [glob_el_id, center_to_find] : serial)	{
		for (const auto& [loc_el_id, local_cent] : local){
			if( center_to_find.DistanceTo(local_cent) <= tol){
				res[glob_el_id] = loc_el_id;
			}
		}
	}
	return res;
}

void ensureElementTypeIsSame(const Mesh& mesh)
{
	const auto type = mesh.GetElementType(0);
	for (auto e = 1; e < mesh.GetNE(); e++){
		if (type != mesh.GetElementType(e)){
			throw std::runtime_error("Parallelization only works if all mesh elements are of the same type");
		}	
	}
}

Model::Model(Mesh& mesh, const GeomTagToMaterialInfo& matInfo, const GeomTagToBoundaryInfo& bdrInfo, int* partitioning)
{

	serialMesh_ = Mesh(mesh);
	ensureElementTypeIsSame(mesh);

	pmesh_ = ParMesh(MPI_COMM_WORLD, serialMesh_, partitioning);

	if (matInfo.gt2m.size() == 0) {
		attToMatMap_.emplace(1, Material(1.0, 1.0, 0.0));
	}
	else {
		attToMatMap_ = matInfo.gt2m;
	}

	if (matInfo.gt2bm.size() != 0){
		attToBdrMatMap_ = matInfo.gt2bm;
	}

	Array<int> f2bdr;
	if (partitioning != nullptr){
		f2bdr = pmesh_.GetFaceToBdrElMap();
	}
	else{
		f2bdr = serialMesh_.GetFaceToBdrElMap();
	}

	for (auto i = bdrInfo.gt2b.begin(); i != bdrInfo.gt2b.end(); i++) {
		faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
	}
	attToBdrMap_ = bdrInfo.gt2b;
	
	if (bdrInfo.gt2ib.size() != 0)
	{
		for (auto i = bdrInfo.gt2ib.begin(); i != bdrInfo.gt2ib.end(); i++){
			if (i->second == BdrCond::PEC || i->second == BdrCond::PMC || i->second == BdrCond::SMA || i->second == BdrCond::SGBC) {
				faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
				attToIntBdrMap_.insert(std::make_pair( i->first, i->second ));
			}
		}
	}
	attToIntBdrMap_ = bdrInfo.gt2ib;

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
	if (sgbc_Marker_.Size() != 0) {
		bdrToMarkerMap_.insert(std::make_pair(BdrCond::SGBC, sgbc_Marker_));
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
	if (intsgbc_Marker_.Size() != 0) {
		intBdrToMarkerMap_.insert(std::make_pair(BdrCond::SGBC, intsgbc_Marker_));
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
			initMarker(getMarker(bdr, isInterior), pmesh_.bdr_attributes.Max());
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
	case BdrCond::SGBC:
		switch (isInterior) {
			case true:
				return intsgbc_Marker_;
			case false:
				return sgbc_Marker_;
		}
		break;
	default:
		throw std::runtime_error("Wrong BdrCond in getMarkerForBdrCond.");
	}
}

}
