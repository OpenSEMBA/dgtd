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
 * Communications.h
 *
 *  Created on: Apr 17, 2013
 *      Author: luis
 */

#ifndef COMMUNICATIONS_H_
#define COMMUNICATIONS_H_

#include <iostream>
#include <assert.h>
#include <vector>

using namespace std;
using namespace SEMBA;
using namespace Geometry;
using namespace Math;

#include "math/Field.h"
#include "Volume.h"
#include "core/Ordering.h"

namespace SEMBA {
namespace Cudg3d {
namespace Communications {

class Comm : public Cell::Ordering {
public:
    virtual ~Comm();
    virtual Int getNumberOfTasks() const = 0;
    virtual Int getTask() const = 0;
    virtual bool isMaster() const = 0;
    virtual Int getNumOfTasksOnThisHost() const = 0;
    virtual size_t getLocalOffset() const = 0;
    virtual size_t getLocalSize() const = 0;
    virtual void setPartitionSizes(const vector<vector<ElemId>>& partId) = 0;
    virtual void gatherFieldsMaster(
            FieldR3& elec, FieldR3& magn,
            const FieldR3& localElec, const FieldR3& localMagn) const = 0;
    virtual void gatherFieldsSlave(
            const FieldR3& electric, const FieldR3& magnetic) const = 0;
    virtual void syncNeighbourFields(
            Real* nEx, Real* nEy, Real* nEz,
            Real* nHx, Real* nHy, Real* nHz,
            const Real* Ex, const Real* Ey, const Real* Ez,
            const Real* Hx, const Real* Hy, const Real* Hz) const = 0;
    virtual void initNeighbourFields(const vector<ElemId>& nIds) = 0;
    virtual Real reduceToGlobalMinimum(Real val) const = 0;
    virtual void printInfo() const = 0;
};

}
}
}

#endif /* COMMUNICATIONS_H_ */
