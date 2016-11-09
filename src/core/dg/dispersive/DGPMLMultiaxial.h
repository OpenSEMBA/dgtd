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
// * SolverPMLUniaxial.h
// *
// *  Created on: Aug 2, 2013
// *      Author: luis
// */
//
//#ifndef SOLVERPMLMULTIAXIAL_H_
//#define SOLVERPMLMULTIAXIAL_H_
//
//#include "DGPML.h"
//
//class DGPMLMultiaxial : public DGPML {
//public:
//    DGPMLMultiaxial(
//            const PMVolumePML& mat,
//            const bool useConductivity,
//            const Real conductivity);
//    virtual ~DGPMLMultiaxial();
//    void addRHSToRes(
//            const size_t e1, const size_t e2,
//            const Real rka, const Real dt);
//    void updateWithRes(
//            const size_t e1,
//            const size_t e2,
//            const Real rkb);
//    virtual void computeRHSElectric(
//            FieldR3& rhsE,
//            const FieldR3& E,
//            const size_t e1, const size_t e2) const = 0;
//    virtual void computeRHSMagnetic(
//            FieldR3& rhsH,
//            const FieldR3& H,
//            const size_t e1, const size_t e2) const = 0;
//    virtual void computeRHSElectricPolarizationCurrents(
//            const FieldR3& E,
//            const size_t e1, const size_t e2) = 0;
//    virtual void computeRHSMagneticPolarizationCurrents(
//            const FieldR3& H,
//            const size_t e1, const size_t e2) = 0;
//protected:
//    FieldR3 J, M, resJ, resM, rhsJ, rhsM;
//};
//
//#endif /* SOLVERPMLUNIAXIAL_H_ */
