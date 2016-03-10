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
 * Communications.cpp
 *
 *  Created on: Apr 17, 2013
 *      Author: luis
 */

#include "../../dgtd/core/CommNone.h"

CommNone::CommNone() {
}

CommNone::~CommNone() {
}

Int CommNone::getNumberOfTasks() const {
   return 1;
}

void CommNone::abort() const {

}

UInt CommNone::getLocalOffset() const {
   return 0;
}

UInt CommNone::getLocalSize() const {
   return getGlobalSize();
}

bool CommNone::isMaster() const {
   return true;
}

void CommNone::gatherFieldsMaster(
      FieldR3& electric,
      FieldR3& magnetic,
      const FieldR3& localElectric,
      const FieldR3& localMagnetic) const {
}

void CommNone::gatherFieldsSlave(
      const FieldR3& electric,
      const FieldR3& magnetic) const {
   //
}

void CommNone::setPartitionSizes(
      const vector<vector<ElemId>>& partId) {
   assert(partId.size() == 1);
   setGlobalSize(partId[0].size());
   setLocalSizeAndOffset(getGlobalSize(), 0);
}

void CommNone::syncNeighbourFields(Real* nEx, Real* nEy, Real* nEz,
      Real* nHx, Real* nHy, Real* nHz, const Real* Ex, const Real* Ey,
      const Real* Ez, const Real* Hx, const Real* Hy,
      const Real* Hz) const {

}

void CommNone::initNeighbourFields(const vector<ElemId>& nIds) {

}

Int CommNone::getTask() const {
   return 0;
}

Real CommNone::reduceToGlobalMinimum(Real val) const {
   return val;
}

Int CommNone::getNumOfTasksOnThisHost() const {
   return 1;
}

void CommNone::printInfo() const {
   cout << "---- CommNone Info ----" <<endl;
   cout << "No multiprocessor comunications being used." << endl;
}
