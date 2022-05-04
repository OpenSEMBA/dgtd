#pragma once

#include "mfem.hpp"
#include "Material.h"
#include <map>

using namespace mfem;

namespace maxwell {

	using attribute = std::size_t;

class Model {
public:

	Model(Mesh& mesh, std::map<attribute, Material>& matMap);
	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; };
	std::map<attribute, Material> getMaterialMap() const { return attToMatMap_; }

private:

	Mesh mesh_;
	std::map<attribute, Material> attToMatMap_;

};

}