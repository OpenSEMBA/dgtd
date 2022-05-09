#pragma once

#include "mfem.hpp"
#include "Material.h"
#include <map>

using namespace mfem;

namespace maxwell {

	using attribute = std::size_t;

class Model {
public:

	Model(Mesh& mesh, std::vector<std::pair<attribute, Material>>& attToMatVec);
	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; };
	const std::vector<std::pair<attribute, Material>>& getAttToMatVec() const { return attToMatVec_; }

private:

	Mesh mesh_;
	std::vector<std::pair<attribute, Material>> attToMatVec_;

};

}