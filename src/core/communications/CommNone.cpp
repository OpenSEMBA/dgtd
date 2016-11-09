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
// * Communications.cpp
// *
// *  Created on: Apr 17, 2013
// *      Author: luis
// */
//
//#include "CommNone.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//namespace Communications {
//
//CommNone::CommNone() {
//}
//
//CommNone::~CommNone() {
//}
//
//Math::Int CommNone::getNumberOfTasks() const {
//   return 1;
//}
//
//void CommNone::abort() const {
//
//}
//
//size_t CommNone::getLocalOffset() const {
//   return 0;
//}
//
//size_t CommNone::getLocalSize() const {
//   return getGlobalSize();
//}
//
//bool CommNone::isMaster() const {
//   return true;
//}
//
//void CommNone::gatherFieldsMaster(
//      Math::FieldR3& electric,
//      Math::FieldR3& magnetic,
//      const Math::FieldR3& localElectric,
//      const Math::FieldR3& localMagnetic) const {
//}
//
//void CommNone::gatherFieldsSlave(
//      const Math::FieldR3& electric,
//      const Math::FieldR3& magnetic) const {
//   //
//}
//
//void CommNone::setPartitionSizes(
//      const vector<vector<Geometry::ElemId>>& partId) {
//   assert(partId.size() == 1);
//   setGlobalSize(partId[0].size());
//   setLocalSizeAndOffset(getGlobalSize(), 0);
//}
//
//void CommNone::syncNeighbourFields(Math::Real* nEx, Math::Real* nEy, Math::Real* nEz,
//      Math::Real* nHx, Math::Real* nHy, Math::Real* nHz, const Math::Real* Ex, const Math::Real* Ey,
//      const Math::Real* Ez, const Math::Real* Hx, const Math::Real* Hy,
//      const Math::Real* Hz) const {
//
//}
//
//void CommNone::initNeighbourFields(const vector<ElemId>& nIds) {
//
//}
//
//Math::Int CommNone::getTask() const {
//   return 0;
//}
//
//Math::Real CommNone::reduceToGlobalMinimum(Math::Real val) const {
//   return val;
//}
//
//Math::Int CommNone::getNumOfTasksOnThisHost() const {
//   return 1;
//}
//
//void CommNone::printInfo() const {
//   cout << "---- CommNone Info ----" <<endl;
//   cout << "No multiprocessor comunications being used." << endl;
//}
//
//}
//}
//}
