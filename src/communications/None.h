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

#ifndef COMMNONE_H_
#define COMMNONE_H_

#include "stdlib.h"
#include <iostream>
#include <vector>
#include <assert.h>

using namespace std;

#include "Communications.h"

namespace SEMBA {
namespace Cudg3d {
namespace Communications {

class None : public Communications {
public:
//    None();
//    virtual ~None();
//    Math::Int getNumberOfTasks() const;
//    void abort() const;
//    bool isMaster() const;
//    Math::Int getTask() const;
//    size_t getLocalOffset() const;
//    Math::Int getNumOfTasksOnThisHost() const;
//    size_t getLocalSize() const;
//    void gatherFieldsMaster(
//            Math::FieldR3& electric,
//            Math::FieldR3& magnetic,
//            const Math::FieldR3& localElectric,
//            const Math::FieldR3& localMagnetic) const;
//    void gatherFieldsSlave(
//            const Math::FieldR3& electric,
//            const Math::FieldR3& magnetic) const;
//    void setPartitionSizes(
//            const vector<vector<ElemId>>& partId);
//    void syncNeighbourFields(
//            Math::Real* nEx, Math::Real* nEy, Math::Real* nEz,
//            Math::Real* nHx, Math::Real* nHy, Math::Real* nHz,
//            const Math::Real* Ex, const Math::Real* Ey, const Math::Real* Ez,
//            const Math::Real* Hx, const Math::Real* Hy, const Math::Real* Hz) const;
//    Math::Real reduceToGlobalMinimum(Math::Real val) const;
//    void initNeighbourFields(const vector<ElemId>& nIds);
//    void printInfo() const;
};

}
}
}

#endif
