#pragma once

#include <mfem.hpp>
#include <map>

#include "Material.h"
#include "Types.h"

namespace maxwell {

using Attribute = int;
using AttributeToMaterial = std::map<Attribute, Material>;
using AttributeToBoundary = std::map<Attribute, BdrCond>;
using AttributeToInteriorBoundary = std::map<Attribute, BdrCond>;

using BoundaryMarker = mfem::Array<int>;
using BoundaryToMarker = std::multimap<BdrCond, BoundaryMarker>;
using InteriorBoundaryToMarker = std::multimap<BdrCond, BoundaryMarker>;

using InteriorBoundaryMarker = mfem::Array<int>;
using InteriorBoundaryToMarker = std::multimap<BdrCond, InteriorBoundaryMarker>;

class Model {
public:
	using Mesh = mfem::Mesh;

	Model(
		Mesh&, 
		const AttributeToMaterial& = AttributeToMaterial{},
		const AttributeToBoundary& = AttributeToBoundary{},
		const AttributeToInteriorBoundary & = AttributeToInteriorBoundary{}
	);

	Mesh& getMesh() { return mesh_; };
	
	BoundaryToMarker& getBoundaryToMarker() { return bdrToMarkerMap_; }
	InteriorBoundaryToMarker& getInteriorBoundaryToMarker() { return intBdrToMarkerMap_; }

	mfem::Vector buildPiecewiseArgVector(const FieldType& f) const;

private:
	Mesh mesh_;
	
	AttributeToMaterial attToMatMap_;
	AttributeToBoundary attToBdrMap_;
	AttributeToInteriorBoundary attToIntBdrMap_;
	BoundaryToMarker bdrToMarkerMap_;
	InteriorBoundaryToMarker intBdrToMarkerMap_;
};

}