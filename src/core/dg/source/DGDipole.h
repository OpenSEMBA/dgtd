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
// * SolverDipole.h
// *
// *  Created on: Oct 3, 2012
// *      Author: luis
// */
//
//#ifndef SOLVERDIPOLE_H_
//#define SOLVERDIPOLE_H_
//
//#include "../../dg/sources/DGSource.h"
//#include "math/SphericalVector.h"
//#include "sources/Dipole.h"
//
//using namespace std;
//
////#define SOLVERDIPOLE_DO_NOT_USE_GAUSSIAN_DERIVATIVE
//#define SOLVERDIPOLE_CREATE_GRAPH_WITH_EXCITATION
//
//class DGDipole : public DGSource, public Dipole {
//public:
//    DGDipole(
//            const Dipole& dip,
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells,
//            FieldR3& dE,
//            FieldR3& dH,
//            const Int vmapM[faces][nfp]);
//    virtual ~DGDipole();
//    void computeExcitation(
//            const Real intTime,
//            const Real minDT);
//    CVecR3 getMagnitude(const Real time) const;
//    void printInfo() const;
//private:
//#ifdef SOLVERDIPOLE_DO_NOT_USE_GAUSSIAN_DERIVATIVE
//    Real *intT, *intS;
//#endif
//    SphericalVector *tPos, *sPos;
//    void computeExcitationField(
//            FieldR3& EInc,
//            FieldR3& HInc,
//            const SphericalVector* vPos,
//            const size_t nE,
//            const Real time) const;
//};
//#endif /* SOLVERDIPOLE_H_ */
