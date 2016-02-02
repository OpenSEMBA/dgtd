// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
// File: Mesh.h
#ifndef MeshVolume_H_
#define MeshVolume_H_

#ifdef USE_METIS
    #include <metis.h>
    #if METIS_VER_MAJOR < 5
        #error "Mesh partitioning requires METIS version 5+"
    #endif
    #define MESH_ALLOW_PARTITIONING
#endif
using namespace std;

#include "geometry/mesh/Unstructured.h"

namespace SEMBA {
namespace Geometry {
namespace Mesh {

class MeshVolume : public Unstructured {
public:
    MeshVolume();
    MeshVolume(const MeshVolume& meshVol);
    MeshVolume(const Unstructured& uns);
    virtual ~MeshVolume();
    MeshVolume& operator=(const MeshVolume& param);
    vector<vector<ElementId>> getPartitionsIds(
            const UInt nDivisions,
            const vector<pair<ElementId,Int>> idWeights = vector<pair<ElementId,Int>>(),
            const Real* taskPower = NULL) const;
};

} /* namespace Mesh */
} /* namespace Geometry */
} /* namespace SEMBA */

#endif
