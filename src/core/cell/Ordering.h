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
// * Ordering.h
// *
// *  Created on: Jun 12, 2013
// *      Author: luis
// */
//
//#ifndef ORDERING_H_
//#define ORDERING_H_
//
//#include "stdlib.h"
//#include <assert.h>
//#include <iostream>
//#include <vector>
//
//using namespace std;
//
//#include "geometry/element/Element.h"
//#include "math/matrix/Dynamic.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//namespace Cell {
//
//class Ordering {
//public:
//    size_t getGlobalSize() const;
//    size_t getGlobalRelPosOfId(const Geometry::ElemId id) const;
//    size_t getRelPosOfId(const Geometry::ElemId id) const;
//    Geometry::ElemId getIdOfGlobalRelPos(const size_t rp) const;
//    Geometry::ElemId getIdOfRelPos(const size_t rp) const;
//    size_t getLocalSize() const;
//    bool checkRelPosOfId() const;
//    bool isLocalId(const Geometry::ElemId id) const;
//protected:
//    Ordering();
//    virtual ~Ordering();
//    void setGlobalSize(const size_t globalSize_);
//    void setLocalSizeAndOffset(
//            const size_t localSize,
//            const size_t localOffset);
//    void buildRelPosOfIds(
//            const Math::Matrix::Dynamic<size_t>& list);
//    void printOrderingInfo() const;
//private:
//    static size_t globalSize;
//    static size_t localSize;
//    static size_t localOffset;
//    static Geometry::ElemId offsetId;
//    static Geometry::ElemId* idOfRelPos;
//    static size_t* relPosOfId;
//    bool checkLocalIds(
//            const vector<vector<Geometry::ElemId> >& partIds,
//            const size_t task);
//};
//
//}
//}
//}
//
//#endif /* ORDERING_H_ */
