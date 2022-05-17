#pragma once

#include "mfem.hpp"
#include "Material.h"
#include <map>

using namespace mfem;

namespace maxwell {

using Attribute = std::size_t;
using AttributeToMaterial = std::map<Attribute, Material>;

class Model {
public:

	Model(Mesh& mesh, const AttributeToMaterial& attToMatVec);
	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; };
	const AttributeToMaterial& getAttToMat() const { return attToMatVec_; }

private:
	Mesh mesh_;
	AttributeToMaterial attToMatVec_;

};

}