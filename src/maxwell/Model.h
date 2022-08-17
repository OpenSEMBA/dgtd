#pragma once

#include <mfem.hpp>
#include <map>

#include "Material.h"
#include "Types.h"

namespace maxwell {

using Attribute = std::size_t;
using BdrIdx = std::size_t;
using AttributeToMaterial = std::map<Attribute, Material>;
using AttributeToBoundary = std::map<Attribute, BdrCond>;

class Model {
public:
	using Mesh = mfem::Mesh;

	Model(Mesh&, const AttributeToMaterial&, const AttributeToBoundary&);
	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; };
	const AttributeToMaterial& getAttToMat() const { return attToMatMap_; }
	const AttributeToBoundary& getAttToBdr() const { return attToBdrMap_; }

private:
	Mesh mesh_;
	AttributeToMaterial attToMatMap_;
	AttributeToBoundary attToBdrMap_;
};

}