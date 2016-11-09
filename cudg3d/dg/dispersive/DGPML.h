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
// * SolverPML.h
// *
// *  Created on: Sep 11, 2012
// *      Author: luis
// */
//
//#ifndef SOLVERPML_H_
//#define SOLVERPML_H_
//
//#include "DGDispersive.h"
//#include "physicalModel/PMVolumePML.h"
//
//class DGPML : public DGDispersive {
//public:
//    DGPML(const PMVolumePML& mat);
//    virtual ~DGPML();
//    void addJumps(
//            FieldR3& dE, FieldR3& dH,
//            FieldR3& E, FieldR3& H,
//            const size_t e1, const size_t e2);
//protected:
//    size_t dof;
//    size_t nElem;
//    size_t *elem;
//    bool useConstantConductivity;
//    static constexpr Real sigDefault = 10e9;
//    Real sig;
//    static const size_t N = ORDER_N;
//    static const size_t np = (N+1) * (N+2) * (N+3) / 6;
//    Real **sig1, **sig2, **sig3;
//    Real **sig11, **sig22, **sig33;
//    Real **sig12, **sig23, **sig31;
//private:
//    void initConductivityMatrices(
//            const PMVolumePML& mat,
//            const CellGroup& cells);
//    void initConductivity(
//            Real **sigma,
//            const size_t ori,
//            const PMVolumePML& mat,
//            const CellGroup& cells);
//};
//
//#endif /* SOLVERPML_H_ */
