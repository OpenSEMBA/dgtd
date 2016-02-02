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
// * SolverPML.cpp
// *
// *  Created on: Sep 11, 2012
// *      Author: luis
// */
//
//#include "DGPML.h"
//
//DGPML::DGPML(PMVolumePML& mat, CellGroup& cells) {
//    // Initializes Element list and dof.
//    vector<ElementId> auxList;
//    auxList.reserve(cells.getLocalSize());
//    for (UInt e = 0; e < cells.getLocalSize(); e++) {
//        const CellTet<ORDER_N>* cell = cells(e);
//        ElementId id = cell->getId();
//        if (cell->material->getId() == mat.getId()) {
//            auxList.push_back(cells.getRelPosOfId(id));
//        }
//        nElem = auxList.size();
//        elem = new UInt[nElem];
//        for (UInt j = 0; j < nElem; j++) {
//            elem[j] = auxList[j];
//        }
//    }
//    dof = np * nElem;
//    // Initializes conductivity matrices.
//    if (!useConstantConductivity) {
//        initConductivityMatrices(mat, cells);
//    }
//}
//
//void DGPML::initConductivity(
//        Real **sigma,
//        const UInt orientation,
//        const PMVolumePML& mat,
//        const CellGroup& cells) {
//    StaMatrix<Real,np,np> auxSig;
//    for (UInt e = 0; e < nElem; e++) {
//        ElementId id = cells.getIdOfRelPos(elem[e]);
//        const CellTet<ORDER_N>* cell;
//        cell = cells.getPtrToCellWithId(id);
//        auxSig =
//                cell->getConductivityWithGeometricProfile(mat, orientation, sig);
//        auxSig.convertToArray(MATRICES_ROW_MAJOR, sigma[e]);
//    }
//}
//
//void DGPML::initConductivityMatrices(
//        const PMVolumePML& mat,
//        const CellGroup& cells) {
//    assert(elem != NULL);
//    if (mat.isUniaxial()) {
//        sig1 = new Real*[nElem];
//        sig11 = new Real*[nElem];
//        for (UInt e = 0; e < nElem; e++) {
//            sig1[e] = new Real[np*np];
//            sig11[e] = new Real[np*np];
//        }
//        initConductivity(sig1, 1, mat, cells);
//        initConductivity(sig11, 11, mat, cells);
//    } else if (mat.isBiaxial()) {
//        sig1 = new Real*[nElem];
//        sig2 = new Real*[nElem];
//        sig11 = new Real*[nElem];
//        sig22 = new Real*[nElem];
//        sig12= new Real*[nElem];
//        for (UInt e = 0; e < nElem; e++) {
//            sig1[e] = new Real[np*np];
//            sig2[e] = new Real[np*np];
//            sig11[e] = new Real[np*np];
//            sig22[e] = new Real[np*np];
//            sig12[e] = new Real[np*np];
//        }
//        initConductivity(sig1, 1, mat, cells);
//        initConductivity(sig2, 2, mat, cells);
//        initConductivity(sig11, 11, mat, cells);
//        initConductivity(sig22, 22, mat, cells);
//        initConductivity(sig12, 12, mat, cells);
//    } else {
//        assert(mat.isTriaxial());
//        sig1 = new Real*[nElem];
//        sig2 = new Real*[nElem];
//        sig3 = new Real*[nElem];
//        sig11 = new Real*[nElem];
//        sig22 = new Real*[nElem];
//        sig33 = new Real*[nElem];
//        sig12= new Real*[nElem];
//        sig23= new Real*[nElem];
//        sig31= new Real*[nElem];
//        for (UInt e = 0; e < nElem; e++) {
//            sig1[e] = new Real[np*np];
//            sig2[e] = new Real[np*np];
//            sig3[e] = new Real[np*np];
//            sig11[e] = new Real[np*np];
//            sig22[e] = new Real[np*np];
//            sig33[e] = new Real[np*np];
//            sig12[e] = new Real[np*np];
//            sig23[e] = new Real[np*np];
//            sig31[e] = new Real[np*np];
//        }
//        initConductivity(sig1, 1, mat, cells);
//        initConductivity(sig2, 2, mat, cells);
//        initConductivity(sig3, 3, mat, cells);
//        initConductivity(sig11, 11, mat, cells);
//        initConductivity(sig22, 22, mat, cells);
//        initConductivity(sig33, 33, mat, cells);
//        initConductivity(sig12, 12, mat, cells);
//        initConductivity(sig23, 23, mat, cells);
//        initConductivity(sig31, 31, mat, cells);
//    }
//}
//
//void DGPML::addJumps(
//        FieldR3& dE, FieldR3& dH,
//        FieldR3& E, FieldR3& H,
//        const UInt e1, const UInt e2) {
//}
