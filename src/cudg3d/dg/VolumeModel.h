#pragma once 

#include "model/Model.h"
#include "geometry/element/Tetrahedron.h"

namespace SEMBA::dgtd::dg {

class VolumeModel {
public:
    VolumeModel(const Model::UnstructuredModel&);

    const Geometry::Mesh::Unstructured& getMesh() const 
    { 
        return mesh_; 
    }

    std::size_t numberOfVolumeElements() const 
    {
        return mesh_.elems().sizeOf<Geometry::Tet>();
    }

private:
    Geometry::Mesh::Unstructured mesh_;
    PMGroup physicalModels_;
};

} 

