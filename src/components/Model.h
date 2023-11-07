#pragma once

#include "Types.h"
#include "Material.h"

#include <mfem.hpp>

#include <map>

namespace maxwell {

using FaceId = int;
using GeomTag = int;
using GeomTagToMaterial = std::map<GeomTag, Material>;
using GeomTagToBoundary = std::map<GeomTag, BdrCond>;
using GeomTagToInteriorConditions = std::map<GeomTag, BdrCond>;
using GeomTagToInteriorSource = std::map<GeomTag, BdrCond>;
using FaceToGeomTag = std::map<FaceId, GeomTag>;
using GeomTagToBdrCond = std::map<GeomTag, BdrCond>;

using BoundaryMarker = mfem::Array<int>;
using InteriorBoundaryMarker = BoundaryMarker;
using BoundaryToMarker = std::map<BdrCond, BoundaryMarker>;
using InteriorBoundaryCondToMarker = std::map<BdrCond, BoundaryMarker>;
using InteriorSourceToMarker = std::map<BdrCond, BoundaryMarker>;

class Model {
public:
	using Mesh = mfem::Mesh;
	Model() = default;
	Model(
		Mesh&, 
		const GeomTagToMaterial& = GeomTagToMaterial{},
		const GeomTagToBoundary& = GeomTagToBoundary{},
		const GeomTagToInteriorConditions& = GeomTagToInteriorConditions{}
	);

	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; }
	
	BoundaryMarker& getMarker(const BdrCond&, bool isInterior);
	
	BoundaryToMarker& getBoundaryToMarker() { return bdrToMarkerMap_; }
	const BoundaryToMarker& getBoundaryToMarker() const { return bdrToMarkerMap_; }
	InteriorBoundaryCondToMarker& getInteriorBoundaryToMarker() { return intBdrToMarkerMap_; }
	InteriorSourceToMarker& getInteriorSourceToMarker() { return intSrcToMarkerMap_; }
	const FaceToGeomTag& getFaceToGeometryTag() { return faceToGeomTag_; }

	mfem::Vector buildPiecewiseArgVector(const FieldType& f) const;

	std::size_t numberOfMaterials() const;
	std::size_t numberOfBoundaryMaterials() const;

private:

	Mesh mesh_;
	
	GeomTagToMaterial attToMatMap_;
	GeomTagToBoundary attToBdrMap_;
	GeomTagToInteriorConditions attToIntBdrMap_;
	GeomTagToInteriorSource attToIntSrcMap_;
	BoundaryToMarker bdrToMarkerMap_;
	InteriorBoundaryCondToMarker intBdrToMarkerMap_;
	InteriorSourceToMarker intSrcToMarkerMap_;
	FaceToGeomTag faceToGeomTag_;

	BoundaryMarker pecMarker_;
	BoundaryMarker pmcMarker_;
	BoundaryMarker smaMarker_;

	BoundaryMarker intpecMarker_;
	BoundaryMarker intpmcMarker_;
	BoundaryMarker intsmaMarker_;

	BoundaryMarker tfsfMarker_;

	void assembleGeomTagToTypeMap(
		std::map<GeomTag, BdrCond>& attToCond, 
		bool isInterior);

	void assembleBdrToMarkerMaps();

	void assembleIfNotEmpty(const BoundaryMarker&, bool isInterior);
};

}