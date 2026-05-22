#include "Model.h"

#include <limits>

#include <unordered_set>

namespace maxwell {

using namespace mfem;

namespace {

std::unordered_set<int> buildPMLTagSet(const GeomTagToPMLRegion& pml_regions)
{
	std::unordered_set<int> tags;
	for (const auto& [tag, _] : pml_regions) {
		tags.insert(tag);
	}
	return tags;
}

std::pair<mfem::Vector, mfem::Vector> computeBoundingBoxForElements(
	const mfem::Mesh& mesh,
	const std::function<bool(int)>& include_attribute)
{
	const int dim = mesh.Dimension();
	mfem::Vector min_corner(dim);
	mfem::Vector max_corner(dim);
	for (int d = 0; d < dim; ++d) {
		min_corner[d] = std::numeric_limits<double>::infinity();
		max_corner[d] = -std::numeric_limits<double>::infinity();
	}

	bool found = false;
	for (int e = 0; e < mesh.GetNE(); ++e) {
		if (!include_attribute(mesh.GetAttribute(e))) {
			continue;
		}

		found = true;
		mfem::Array<int> vertices;
		mesh.GetElementVertices(e, vertices);
		for (int i = 0; i < vertices.Size(); ++i) {
			const double* coords = mesh.GetVertex(vertices[i]);
			for (int d = 0; d < dim; ++d) {
				min_corner[d] = std::min(min_corner[d], coords[d]);
				max_corner[d] = std::max(max_corner[d], coords[d]);
			}
		}
	}

	if (!found) {
		throw std::runtime_error("Could not infer a PML box because no matching elements were found.");
	}

	return std::make_pair(min_corner, max_corner);
}

}


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

Model::Model(Mesh& mesh, const GeomTagToMaterialInfo& matInfo, const GeomTagToBoundaryInfo& bdrInfo, int* partitioning, MPI_Comm comm)
{

	serialMesh_ = Mesh(mesh);
	ensureElementTypeIsSame(mesh);

	pmesh_ = ParMesh(comm, serialMesh_, partitioning);

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
	const std::pair<BdrCond, const mfem::Array<int>&> bdrEntries[] = {
		{BdrCond::PEC, pecMarker_}, {BdrCond::PMC, pmcMarker_},
		{BdrCond::SMA, smaMarker_}, {BdrCond::SGBC, sgbc_Marker_}
	};
	for (const auto& [cond, marker] : bdrEntries) {
		if (marker.Size() != 0) bdrToMarkerMap_.insert({cond, marker});
	}

	const std::pair<BdrCond, const mfem::Array<int>&> intBdrEntries[] = {
		{BdrCond::PEC, intpecMarker_}, {BdrCond::PMC, intpmcMarker_},
		{BdrCond::SMA, intsmaMarker_}, {BdrCond::SGBC, intsgbc_Marker_}
	};
	for (const auto& [cond, marker] : intBdrEntries) {
		if (marker.Size() != 0) intBdrToMarkerMap_.insert({cond, marker});
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

mfem::Array<int> Model::buildPMLVolumeMarker() const
{
	mfem::Array<int> marker;
	if (pml_regions_.empty()) {
		return marker;
	}

	marker.SetSize(pmesh_.attributes.Max());
	marker = 0;
	for (const auto& [tag, _] : pml_regions_) {
		if (tag <= 0 || tag > marker.Size()) {
			throw std::runtime_error("PML geometry tag " + std::to_string(tag) + " is outside the mesh attribute range.");
		}
		marker[tag - 1] = 1;
	}

	return marker;
}

PMLBoxGeometry Model::inferPMLBoxGeometry() const
{
	if (pml_regions_.empty()) {
		throw std::runtime_error("Cannot infer PML box geometry without any PML regions.");
	}

	const auto pml_tags = buildPMLTagSet(pml_regions_);
	const auto outer_bbox = computeBoundingBoxForElements(serialMesh_, [](int) { return true; });
	const auto inner_bbox = computeBoundingBoxForElements(
		serialMesh_,
		[&pml_tags](int attr) {
			return pml_tags.find(attr) == pml_tags.end();
		});

	PMLBoxGeometry geom;
	const int dim = serialMesh_.Dimension();
	geom.outer_min = outer_bbox.first;
	geom.outer_max = outer_bbox.second;
	geom.inner_min = inner_bbox.first;
	geom.inner_max = inner_bbox.second;
	geom.thickness_minus.SetSize(dim);
	geom.thickness_plus.SetSize(dim);
	bool has_active_axis = false;

	for (int d = 0; d < dim; ++d) {
		const double span = std::max(geom.outer_max[d] - geom.outer_min[d], 1.0);
		const double tol = 1e-8 * span;
		const double minus = geom.inner_min[d] - geom.outer_min[d];
		const double plus = geom.outer_max[d] - geom.inner_max[d];

		if (minus < -tol || plus < -tol) {
			throw std::runtime_error("PML shell geometry is invalid: inner box extends outside the full mesh bounding box.");
		}

		geom.thickness_minus[d] = std::max(0.0, minus);
		geom.thickness_plus[d] = std::max(0.0, plus);
		geom.active_axes[d] = geom.thickness_minus[d] > tol || geom.thickness_plus[d] > tol;
		has_active_axis = has_active_axis || geom.active_axes[d];

		if (geom.active_axes[d] && (geom.thickness_minus[d] <= tol || geom.thickness_plus[d] <= tol)) {
			throw std::runtime_error("PML shell geometry is invalid: active PML axes must wrap both sides of the inner box.");
		}
	}

	if (!has_active_axis) {
		throw std::runtime_error("PML shell geometry is invalid: the tagged region does not define an outer shell around the inner box.");
	}

	return geom;
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
