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

#ifndef COMMNONE_H_
#define COMMNONE_H_

#include "stdlib.h"
#include <iostream>
#include <vector>
#include <assert.h>

#include "../../dgtd/core/Comm.h"

using namespace std;

class CommNone : public Comm {
public:
    CommNone();
    virtual ~CommNone();
    Int getNumberOfTasks() const;
    void abort() const;
    bool isMaster() const;
    Int getTask() const;
    UInt getLocalOffset() const;
    Int getNumOfTasksOnThisHost() const;
    UInt getLocalSize() const;
    void gatherFieldsMaster(
            FieldR3& electric,
            FieldR3& magnetic,
            const FieldR3& localElectric,
            const FieldR3& localMagnetic) const;
    void gatherFieldsSlave(
            const FieldR3& electric,
            const FieldR3& magnetic) const;
    void setPartitionSizes(
            const vector<vector<ElemId>>& partId);
    void syncNeighbourFields(
            Real* nEx, Real* nEy, Real* nEz,
            Real* nHx, Real* nHy, Real* nHz,
            const Real* Ex, const Real* Ey, const Real* Ez,
            const Real* Hx, const Real* Hy, const Real* Hz) const;
    Real reduceToGlobalMinimum(Real val) const;
    void initNeighbourFields(const vector<ElemId>& nIds);
    void printInfo() const;
};

#endif
