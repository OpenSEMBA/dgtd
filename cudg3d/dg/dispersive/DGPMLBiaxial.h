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
// * DGPMLBiaxial.h
// *
// *  Created on: Jun 22, 2015
// *      Author: luis
// */
//
//#ifndef SRC_SOLVER_DGTD_DG_DISPERSIVES_DGPMLBIAXIAL_H_
//#define SRC_SOLVER_DGTD_DG_DISPERSIVES_DGPMLBIAXIAL_H_
//
//#include "DGPMLMultiaxial.h"
//
//template<Int D>
//class DGPMLBiaxial: public DGPMLMultiaxial {
//public:
//    DGPMLBiaxial(
//            const PMVolumePML& mat,
//            const CellGroup& cells,
//            const bool useConductivity,
//            const Real conductivity);
//    virtual ~DGPMLBiaxial();
//    void computeRHSElectric(
//            FieldR3& rhs,
//            const FieldR3& f,
//            const size_t e1, const size_t e2) const;
//    void computeRHSMagnetic(
//            FieldR3& rhs,
//            const FieldR3& f,
//            const size_t e1, const size_t e2) const;
//    void computeRHSElectricPolarizationCurrents(
//            const FieldR3& f,
//            const size_t e1, const size_t e2);
//    void computeRHSMagneticPolarizationCurrents(
//            const FieldR3& f,
//            const size_t e1, const size_t e2);
//};
//
//#include "DGPMLBiaxial.hpp"
//
//typedef DGPMLBiaxial<x> DGPMLxy;
//typedef DGPMLBiaxial<y> DGPMLyz;
//typedef DGPMLBiaxial<z> DGPMLzx;
//
//#endif /* SRC_SOLVER_DGTD_DG_DISPERSIVES_DGPMLBIAXIAL_H_ */
