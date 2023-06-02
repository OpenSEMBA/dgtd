#pragma once

#include "Types.h"
#include "Material.h"

#include <mfem.hpp>

#include <map>

namespace maxwell {

using Attribute = int;
using AttributeToMaterial = std::map<Attribute, Material>;
using AttributeToBoundary = std::map<Attribute, BdrCond>;
using AttributeToInteriorConditions = std::map<Attribute, BdrCond>;
using AttributeToInteriorSource = std::map<Attribute, BdrCond>;

using BoundaryMarker = mfem::Array<int>;
using BoundaryToMarker = std::multimap<BdrCond, BoundaryMarker>;
using InteriorBoundaryCondToMarker = std::multimap<BdrCond, BoundaryMarker>;
using InteriorSourceToMarker = std::multimap<BdrCond, BoundaryMarker>;

using InteriorBoundaryMarker = mfem::Array<int>;
using InteriorBoundaryCondToMarker = std::multimap<BdrCond, InteriorBoundaryMarker>;

class Model {
public:
	using Mesh = mfem::Mesh;
	Model() = default;
	Model(
		Mesh&, 
		const AttributeToMaterial& = AttributeToMaterial{},
		const AttributeToBoundary& = AttributeToBoundary{},
		const AttributeToInteriorConditions & = AttributeToInteriorConditions{}
	);

	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; }
	
	BoundaryToMarker& getBoundaryToMarker() { return bdrToMarkerMap_; }
	const BoundaryToMarker& getBoundaryToMarker() const { return bdrToMarkerMap_; }
	InteriorBoundaryCondToMarker& getInteriorBoundaryToMarker() { return intBdrToMarkerMap_; }
	InteriorSourceToMarker& getInteriorSourceToMarker() { return intSrcToMarkerMap_; }

	mfem::Vector buildPiecewiseArgVector(const FieldType& f) const;

	std::size_t numberOfMaterials() const;
	std::size_t numberOfBoundaryMaterials() const;
private:
	Mesh mesh_;
	
	AttributeToMaterial attToMatMap_;
	AttributeToBoundary attToBdrMap_;
	AttributeToInteriorConditions attToIntBdrMap_;
	AttributeToInteriorSource attToIntSrcMap_;
	BoundaryToMarker bdrToMarkerMap_;
	InteriorBoundaryCondToMarker intBdrToMarkerMap_;
	InteriorSourceToMarker intSrcToMarkerMap_;

	void assembleAttToTypeMap(
		std::map<Attribute, BdrCond>& attToCond, 
		std::multimap<BdrCond, BoundaryMarker>& attToMarker);
};

}