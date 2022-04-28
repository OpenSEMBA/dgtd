#pragma once

#include "mfem.hpp"
#include "Material.h"
#include <map>

using namespace mfem;

namespace maxwell {

class Model {
public:
	using attribute = std::size_t;
	using attToMaterialMap = std::map<attribute, Material>;

	const Model(Mesh mesh, attToMaterialMap matMap);
	Mesh& getMesh() { return mesh_; };
	attToMaterialMap getMaterialMap() const { return matMap_; }

private:

	Mesh mesh_;
	attToMaterialMap matMap_;


};

}