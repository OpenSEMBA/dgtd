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
#ifndef MESH_VOLUME_H_
#define MESH_VOLUME_H_

#include <metis.h>
#if METIS_VER_MAJOR < 5
#error "Mesh partitioning requires METIS version 5+"
#endif
//#define MESH_ALLOW_PARTITIONING

using namespace std;

#include "geometry/mesh/Unstructured.h"

namespace SEMBA {
namespace Cudg3d {
namespace Mesh {

class Volume : public SEMBA::Geometry::Mesh::Unstructured {
public:
    Volume();
    Volume(const Volume& meshVol);
    Volume(const Unstructured& uns);
    virtual ~Volume();
    Volume& operator=(const Volume& param);
    vector<vector<Geometry::ElemId>> getPartitionsIds(
            const size_t nDivisions,
            const vector<pair<Geometry::ElemId,int>> idWeights =
                    vector<pair<Geometry::ElemId,int>>(),
            const Math::Real* taskPower = NULL) const;

    const Geometry::Graph::Connectivities* getConnectivities() const;
private:
    Geometry::Graph::Connectivities* connectivities_;
};

} /* namespace Mesh */
}
}

#endif
