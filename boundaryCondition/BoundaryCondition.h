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

#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>			// Stream I/O.
#include <cmath>
#include <vector>
#include <utility>

using namespace std;

#include "Data.h"
#include "cell/Group.h"
#include "cell/Cell.h"

namespace SEMBA {
namespace Cudg3d {
namespace BoundaryCondition {

template<class T>
class BoundaryCondition {
public:
    BoundaryCondition(
            T* condition,
            Cell::Face localFace,
            Cell::Face neighFace);
    virtual ~BoundaryCondition();
    bool hasSameBoundary(const BoundaryCondition& other) const;
    virtual BoundaryCondition& operator=(const BoundaryCondition& rhs);

    Cell::Face getLocalFace() const;
    Cell::Face getNeighFace() const;
    const T* getCondition() const;

    void printInfo() const;
private:
    const T* condition_;
    Cell::Face localFace_, neighFace_;
};

#include "BoundaryCondition.hpp"

typedef BoundaryCondition<Source::Base> EMSourceBC;
typedef BoundaryCondition<PhysicalModel::PhysicalModel> PhysicalModelBC;
typedef BoundaryCondition<PhysicalModel::Surface::SIBC> SurfaceImpedanceBC;

}
}
}


#endif /* BOUNDARYCONDITION_H_ */
