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
// * SolverDispersive.h
// *
// *  Created on: Sep 27, 2014
// *      Author: luis
// */
//
//#ifndef SOLVERDISPERSIVE_H_
//#define SOLVERDISPERSIVE_H_
//
//#include <cmath>
//
//using namespace std;
//
//#include "math/Constants.h"
//#include "math/Field.h"
//#include "physicalModel/PhysicalModel.h"
//#include "solver/dgtd/core/CellGroup.h"
//
//class DGDispersive {
//public:
//    static const size_t N = ORDER_N;
//    static const size_t nfp = (N+1) * (N+2) / 2;
//    static const size_t np = (N+1) * (N+2) * (N+3) / 6;
//    static const size_t faces = 4;
//    virtual ~DGDispersive() = 0;
//    virtual void computeRHSElectric(
//            FieldR3& rhsE,
//            const FieldR3& E,
//            const size_t e1, const size_t e2) const = 0;
//    virtual void computeRHSMagnetic(
//            FieldR3& rhsH,
//            const FieldR3& H,
//            const size_t e1, const size_t e2) const= 0;
//    virtual void computeRHSElectricPolarizationCurrents(
//            const FieldR3& E,
//            const size_t e1, const size_t e2) = 0;
//    virtual void computeRHSMagneticPolarizationCurrents(
//            const FieldR3& H,
//            const size_t e1, const size_t e2) = 0;
//    virtual void addRHSToRes(
//            const size_t e1, const size_t e2,
//            const Real rka, const Real dt) = 0;
//    virtual void updateWithRes(
//            const size_t e1,
//            const size_t e2,
//            const Real rkb) = 0;
//    virtual void addJumps(
//            FieldR3& dE, FieldR3& dH,
//            FieldR3& E, FieldR3& H,
//            const size_t e1, const size_t e2) = 0;
//};
//
//#endif /* SOLVERDISPERSIVE_H_ */
