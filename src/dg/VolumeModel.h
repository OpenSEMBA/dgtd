#pragma once 

#include "geometry/mesh/Unstructured.h"

namespace SEMBA::dgtd::dg {

class VolumeModel {
public:
    VolumeModel(const Geometry::Mesh::Unstructured& uns);

private:
    Geometry::Mesh::Unstructured mesh_;
    PMGroup physicalModels_;
};

} 

