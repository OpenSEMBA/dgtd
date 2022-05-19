#pragma once

#include "mfem.hpp"
#include "Material.h"
#include "Types.h"
#include <map>

using namespace mfem;

namespace maxwell {

using Attribute = std::size_t;
using BdrIdx = std::size_t;
using AttributeToMaterial = std::map<Attribute, Material>;
using AttributeToBoundary = std::map<Attribute, BdrCond>;

class Model {
public:

	Model(Mesh& mesh, const AttributeToMaterial& attToMatVec, const AttributeToBoundary& bdrVec);
	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; };
	const AttributeToMaterial& getAttToMat() const { return attToMatVec_; }
	const AttributeToBoundary& getAttToBdr() const { return attToBdrVec_; }

private:
	Mesh mesh_;
	AttributeToMaterial attToMatVec_;
	AttributeToBoundary attToBdrVec_;

};

}