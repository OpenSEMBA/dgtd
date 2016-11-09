//// OpenSEMBA
//// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
////                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
////                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
////                    Daniel Mateos Romero            (damarro@semba.guru)
////
//// This file is part of OpenSEMBA.
////
//// OpenSEMBA is free software: you can redistribute it and/or modify it under
//// the terms of the GNU Lesser General Public License as published by the Free
//// Software Foundation, either version 3 of the License, or (at your option)
//// any later version.
////
//// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
//// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//// details.
////
//// You should have received a copy of the GNU Lesser General Public License
//// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
///*
// * Communications.h
// *
// *  Created on: Apr 17, 2013
// *      Author: luis
// */
//
//#ifndef COMMUNICATIONS_H_
//#define COMMUNICATIONS_H_
//
//#include <iostream>
//#include <assert.h>
//#include <vector>
//
//using namespace std;
//
//#include "math/Field.h"
//#include "mesh/Volume.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//namespace Communications {
//
//class Comm {
//public:
//    virtual ~Comm();
//    virtual Math::Int getNumberOfTasks() const = 0;
//    virtual Math::Int getTask() const = 0;
//    virtual bool isMaster() const = 0;
//    virtual Math::Int getNumOfTasksOnThisHost() const = 0;
//    virtual size_t getLocalOffset() const = 0;
//    virtual size_t getLocalSize() const = 0;
//    virtual void setPartitionSizes(const vector<vector<Geometry::ElemId>>& partId) = 0;
//    virtual void gatherFieldsMaster(
//            Math::FieldR3& elec, Math::FieldR3& magn,
//            const Math::FieldR3& localElec, const Math::FieldR3& localMagn) const = 0;
//    virtual void gatherFieldsSlave(
//            const Math::FieldR3& electric, const Math::FieldR3& magnetic) const = 0;
//    virtual void syncNeighbourFields(
//            Math::Real* nEx, Math::Real* nEy, Math::Real* nEz,
//            Math::Real* nHx, Math::Real* nHy, Math::Real* nHz,
//            const Math::Real* Ex, const Math::Real* Ey, const Math::Real* Ez,
//            const Math::Real* Hx, const Math::Real* Hy, const Math::Real* Hz) const = 0;
//    virtual void initNeighbourFields(const vector<ElemId>& nIds) = 0;
//    virtual Math::Real reduceToGlobalMinimum(Math::Real val) const = 0;
//    virtual void printInfo() const = 0;
//};
//
//}
//}
//}
//
//#endif /* COMMUNICATIONS_H_ */
