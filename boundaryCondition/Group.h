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

/*
 * BCGroup.h
 *
 *  Created on: Jul 8, 2013
 *      Author: luis
 */

#ifndef BCGROUP_H_
#define BCGROUP_H_

using namespace std;

#include "BoundaryCondition.h"
#include "mesh/Volume.h"
#include "source/Group.h"
#include "physicalModel/Group.h"
#include "geometry/graph/Connectivities.h"
#include "group/Cloneable.h"
#include "group/Printable.h"

namespace SEMBA {
namespace Cudg3d {
namespace BoundaryCondition {

class Group : public SEMBA::Group::Group<Base> {
public:
    Group(const Mesh::Volume& mesh,
            const SourceGroup& em,
            const PMGroup& pm);
    void printInfo() const;

private:
    void buildEMSourceBC_(
            const Mesh::Volume& mesh,
            const SourceGroup& em);
    void buildPhysicalModelBC_(
            const Mesh::Volume& mesh,
            const PMGroup& pm);
    void removeOverlapped();
    vector<Base*> removeCommons(
            const vector<Base*>& low,
            const vector<Base*>& high) const;
};

} /* SEMBA */

typedef BoundaryCondition::Group BCGroup;

} /* Cudg3d */
} /* BoundaryCondition */


#endif /* BCGROUP_H_ */
