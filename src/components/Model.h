#pragma once

#include "Types.h"
#include "Material.h"

namespace maxwell {

using namespace mfem;

using LocalElementId = int;
using NeighbourElementId = int;
using GlobalElementId = int;
using Position = Vector;

using GlobalToLocalElMap = std::map<GlobalElementId, LocalElementId>;

using FaceId = int;
using GeomTag = int;
using GeomTagToMaterial = std::map<GeomTag, Material>;
using GeomTagToBoundaryMaterial = std::map<GeomTag, Material>;
using GeomTagToBoundary = std::map<GeomTag, BdrCond>;
using GeomTagToInteriorBoundary = std::map<GeomTag, BdrCond>;
using GeomTagToInteriorSource = std::map<GeomTag, BdrCond>;
using FaceToGeomTag = std::map<FaceId, GeomTag>;
using GeomTagToBdrCond = std::map<GeomTag, BdrCond>;

using BoundaryMarker = mfem::Array<int>;
using InteriorBoundaryMarker = BoundaryMarker;
using BoundaryToMarker = std::map<BdrCond, BoundaryMarker>;
using InteriorBoundaryToMarker = BoundaryToMarker;
using TotalFieldScatteredFieldToMarker = BoundaryToMarker;
using SGBCToMarker = BoundaryToMarker;
using InteriorSourceToMarker = BoundaryToMarker;

struct GeomTagToBoundaryInfo {
	GeomTagToBoundary gt2b;
	GeomTagToInteriorBoundary gt2ib;

	GeomTagToBoundaryInfo() 
	{
		gt2b = GeomTagToBoundary{};
		gt2ib = GeomTagToInteriorBoundary{};
	};

	GeomTagToBoundaryInfo(const GeomTagToBoundary& gttb, const GeomTagToInteriorBoundary& gttib)
	{
		gt2b = gttb;
		gt2ib = gttib;
	}
};

struct GeomTagToMaterialInfo {
	GeomTagToMaterial gt2m;
	GeomTagToBoundaryMaterial gt2bm;

	GeomTagToMaterialInfo()
	{
		gt2m = GeomTagToMaterial{};
		gt2bm = GeomTagToBoundaryMaterial{};
	};	
	
	GeomTagToMaterialInfo(const GeomTagToMaterial& gttm, const GeomTagToBoundaryMaterial& gttbm)
	{
		gt2m = gttm;
		gt2bm = gttbm;
	};
};

struct SGBCBoundaryInfo
{
    BdrCond bdrCond = BdrCond::SMA;
    bool isOn = false;
};

using SGBCBoundaries = std::pair<SGBCBoundaryInfo, SGBCBoundaryInfo>;

struct SGBCLayer {
    Material material;
    double width;
    size_t num_of_segments = 10;
    size_t order = 1;
    double n_skin_depths = 0.0;

    SGBCLayer(Material mat, double w) : material(mat), width(w) {}
};

struct SGBCProperties{

    std::vector<size_t> geom_tags;
    std::vector<SGBCLayer> layers;
	SGBCBoundaries sgbc_bdr_info;
	bool exporter_probe = false;

    SGBCProperties() : layers() {}

    // Convenience: total width across all layers
    double totalWidth() const {
        double w = 0.0;
        for (const auto& l : layers) w += l.width;
        return w;
    }

    // Convenience: total number of segments across all layers
    size_t totalSegments() const {
        size_t n = 0;
        for (const auto& l : layers) n += l.num_of_segments;
        return n;
    }

    // Convenience: maximum order across all layers
    size_t maxOrder() const {
        size_t o = 1;
        for (const auto& l : layers) o = std::max(o, l.order);
        return o;
    }
};

struct PMLRegionProperties {
	bool matches_vacuum = true;
	size_t grading_order = 3;
	double target_reflection = 1e-6;
	double kappa_max = 1.0;
	double alpha_max = 0.0;

	bool operator==(const PMLRegionProperties& other) const
	{
		return matches_vacuum == other.matches_vacuum &&
			   grading_order == other.grading_order &&
			   target_reflection == other.target_reflection &&
			   kappa_max == other.kappa_max &&
			   alpha_max == other.alpha_max;
	}
};

using GeomTagToPMLRegion = std::map<GeomTag, PMLRegionProperties>;

struct PMLBoxGeometry {
	mfem::Vector inner_min;
	mfem::Vector inner_max;
	mfem::Vector outer_min;
	mfem::Vector outer_max;
	mfem::Vector thickness_minus;
	mfem::Vector thickness_plus;
	std::array<bool, 3> active_axes{ false, false, false };
};

std::map<GlobalElementId, Position> buildSerialElem2CenterMap(Mesh&);
std::map<LocalElementId, Position> buildPartitionElem2CenterMap(ParMesh&);
GlobalToLocalElMap buildGlobalToPartitionLocalElementMap(
	const std::map<GlobalElementId, Position>& serial, const std::map<LocalElementId, Position>& local);
void ensureElementTypeIsSame(const Mesh& mesh);

class Model {
public:

	Model() = default;
	Model(
		Mesh&, 
		const GeomTagToMaterialInfo& = GeomTagToMaterialInfo{},
		const GeomTagToBoundaryInfo& = GeomTagToBoundaryInfo{},
		int* partitioning = nullptr,
		MPI_Comm comm = MPI_COMM_WORLD
	);

	ParMesh& getMesh() { return pmesh_; };
	const ParMesh& getConstMesh() const { return pmesh_; }
	Mesh& getSerialMesh() { return serialMesh_; }
	const Mesh& getConstSerialMesh() const { return serialMesh_; }
	
	BoundaryMarker& getMarker(const BdrCond&, bool isInterior);
	
	const GeomTagToMaterial& getGeomTagToMaterial() { return attToMatMap_; }
	const GeomTagToBoundaryMaterial& getGeomTagToBoundaryMaterial() { return attToBdrMatMap_; }
	BoundaryToMarker& getBoundaryToMarker() { return bdrToMarkerMap_; }
	const BoundaryToMarker& getBoundaryToMarker() const { return bdrToMarkerMap_; }
	InteriorBoundaryToMarker& getInteriorBoundaryToMarker() { return intBdrToMarkerMap_; }
	const InteriorBoundaryToMarker& getInteriorBoundaryToMarker() const { return intBdrToMarkerMap_; }
	TotalFieldScatteredFieldToMarker& getTotalFieldScatteredFieldToMarker() { return tfsfToMarkerMap_; }
	SGBCToMarker& getSGBCToMarker() { return SGBCToMarkerMap_; }
	InteriorSourceToMarker& getInteriorSourceToMarker() { return intSrcToMarkerMap_; }
	const FaceToGeomTag& getFaceToGeometryTag() { return faceToGeomTag_; }
	GeomTagToInteriorBoundary& getGeomTagToIntBoundaryCond() { return attToIntBdrMap_; }
	const GeomTagToInteriorBoundary& getGeomTagToIntBoundaryCond() const { return attToIntBdrMap_; }	
	GeomTagToBoundary& getGeomTagToBoundaryCond() { return attToBdrMap_; }
	const GeomTagToBoundary& getGeomTagToBoundaryCond() const { return attToBdrMap_; }

	void setSGBCProperties(const std::vector<SGBCProperties> in) { sgbc_props_ = in; }
	const std::vector<SGBCProperties>& getSGBCProperties() const { return sgbc_props_; }
	void setPMLRegions(const GeomTagToPMLRegion& in) { pml_regions_ = in; }
	const GeomTagToPMLRegion& getPMLRegions() const { return pml_regions_; }
	bool hasPMLVolumes() const { return !pml_regions_.empty(); }
	mfem::Array<int> buildPMLVolumeMarker() const;
	PMLBoxGeometry inferPMLBoxGeometry() const;

	mfem::Vector initialiseGeomTagVector() const;
	mfem::Vector buildEpsMuPiecewiseVector(const FieldType& f) const;
	mfem::Vector buildSigmaPiecewiseVector() const;

	std::size_t numberOfMaterials() const;
	std::size_t numberOfBoundaryMaterials() const;

	std::string meshName_;

private:

	Mesh serialMesh_;
	ParMesh pmesh_;
	
	GeomTagToMaterial attToMatMap_;
	GeomTagToBoundaryMaterial attToBdrMatMap_;
	GeomTagToBoundary attToBdrMap_;
	GeomTagToInteriorBoundary attToIntBdrMap_;
	GeomTagToInteriorSource attToIntSrcMap_;
	BoundaryToMarker bdrToMarkerMap_;
	InteriorBoundaryToMarker intBdrToMarkerMap_;
	TotalFieldScatteredFieldToMarker tfsfToMarkerMap_;
	SGBCToMarker SGBCToMarkerMap_;
	InteriorSourceToMarker intSrcToMarkerMap_;
	FaceToGeomTag faceToGeomTag_;

	BoundaryMarker pecMarker_;
	BoundaryMarker pmcMarker_;
	BoundaryMarker smaMarker_;

	BoundaryMarker intpecMarker_;
	BoundaryMarker intpmcMarker_;
	BoundaryMarker intsmaMarker_;

	BoundaryMarker tfsfMarker_;

	BoundaryMarker sgbc_Marker_;
	BoundaryMarker intsgbc_Marker_;
	std::vector<SGBCProperties> sgbc_props_;
	GeomTagToPMLRegion pml_regions_;

	void assembleGeomTagToTypeMap(
		std::map<GeomTag, BdrCond>& attToCond, 
		bool isInterior);

	void assembleBdrToMarkerMaps();
};

}