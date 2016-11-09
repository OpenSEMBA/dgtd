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
// * SolverDispersive.cpp
// *
// * Created on: Sep 11, 2012
// *   Author: luis
// */
//
//#include "DGDispersiveVolumic.h"
//
//#include "../DG.h"
//
//DGDispersiveVolumic::DGDispersiveVolumic() {
//    nElem = 0;
//    dof = 0;
//    elem = NULL;
//    drudeDof = 0;
//}
//
//DGDispersiveVolumic::DGDispersiveVolumic(const PMVolumeDispersive& mat) :
//                 PMVolumeDispersive(mat) {
//    //   build(cells);
//}
//
//DGDispersiveVolumic::~DGDispersiveVolumic() {
//
//}
//
//void DGDispersiveVolumic::updateWithRes(
//        const size_t e1,
//        const size_t e2,
//        const Real rkb) {
//    size_t i, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,e)
//#endif
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            P.set(x)[i] += resP(x)[i] * rkb;
//            P.set(y)[i] += resP(y)[i] * rkb;
//            P.set(z)[i] += resP(z)[i] * rkb;
//        }
//    }
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,e)
//#endif
//    for (i=0; i < drudeDof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            J.set(x)[i] += resJ(x)[i] * rkb;
//            J.set(y)[i] += resJ(y)[i] * rkb;
//            J.set(z)[i] += resJ(z)[i] * rkb;
//        }
//    }
//}
//
//void DGDispersiveVolumic::addRHSToRes(
//        const size_t e1,
//        const size_t e2,
//        const Real rka,
//        const Real dt) {
//    size_t i, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,e)
//#endif
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            resP.set(x)[i] *= rka;
//            resP.set(y)[i] *= rka;
//            resP.set(z)[i] *= rka;
//            resP.set(x)[i] += rhsP(x)[i] * dt;
//            resP.set(y)[i] += rhsP(y)[i] * dt;
//            resP.set(z)[i] += rhsP(z)[i] * dt;
//        }
//    }
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,e)
//#endif
//    for (i = 0; i < drudeDof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            resJ.set(x)[i] *= rka;
//            resJ.set(y)[i] *= rka;
//            resJ.set(z)[i] *= rka;
//            resJ.set(x)[i] += rhsJ(x)[i] * dt;
//            resJ.set(y)[i] += rhsJ(y)[i] * dt;
//            resJ.set(z)[i] += rhsJ(z)[i] * dt;
//        }
//    }
//}
//
//void DGDispersiveVolumic::computeRHSElectric(
//        FieldR3& rhsE,
//        const FieldR3& E,
//        const size_t e1, const size_t e2) const {
//    size_t i, j, e, n, p;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n,p)
//#endif
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem]; // Element number.
//        if (e1 <= e && e < e2) {
//            p = i / (nElem * np); // Pole number.
//            n = i % np; // Node number.
//            j = e * np + n; // Field coeff pos.
//            rhsE.set(x)[j] -= 2.0 * real(getPole(p) * P(x)[i] + getResidue(p) * E(x)[j]);
//            rhsE.set(y)[j] -= 2.0 * real(getPole(p) * P(y)[i] + getResidue(p) * E(y)[j]);
//            rhsE.set(z)[j] -= 2.0 * real(getPole(p) * P(z)[i] + getResidue(p) * E(z)[j]);
//        }
//    }
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n,p)
//#endif
//    for (i = 0; i < drudeDof; i++) {
//        // Pole number.
//        p = i / (nElem * np);
//        // Element number.
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            // Node number.
//            n = i % np;
//            // Field coefficient position in the general fields vector.
//            j = e * np + n;
//            // Adds polarization current contribution to the RHS.
//            rhsE.set(x)[j] -= 2.0 * real(J(x)[i]);
//            rhsE.set(y)[j] -= 2.0 * real(J(y)[i]);
//            rhsE.set(z)[j] -= 2.0 * real(J(z)[i]);
//        }
//    }
//}
//
//void DGDispersiveVolumic::computeRHSMagnetic(
//        FieldR3& rhsE,
//        const FieldR3& E, const size_t e1, const size_t e2) const {
//}
//
//void DGDispersiveVolumic::computeRHSElectricPolarizationCurrents(
//        const FieldR3& E,
//        const size_t e1, const size_t e2) {
//    size_t i, j, e, n, p;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n,p)
//#endif
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem]; // Element number.
//        if (e1 <= e && e < e2) {
//            p = i / (nElem * np); // Pole number.
//            n = i % np; // Node number.
//            j = e * np + n; // Field coeff pos.
//            rhsP.set(x)[i] = getPole(p) * P(x)[i] + getResidue(p) * E(x)[j];
//            rhsP.set(y)[i] = getPole(p) * P(y)[i] + getResidue(p) * E(y)[j];
//            rhsP.set(z)[i] = getPole(p) * P(z)[i] + getResidue(p) * E(z)[j];
//        }
//    }
//}
//
//
//void DGDispersiveVolumic::computeRHSMagneticPolarizationCurrents(
//        const FieldR3& E, const size_t e1, const size_t e2) {
//}
//
//void DGDispersiveVolumic::addJumps(FieldR3& dE,
//        FieldR3& dH, FieldR3& E, FieldR3& H,
//        const size_t e1, const size_t e2) {
//}
//
//void DGDispersiveVolumic::build() {
////    for (size_t p = 0; p < getPoleNumber(); p++) {
////        getResidue(p) *= Constants::eps0;
////    }
////    for (size_t p = 0; p < getDrudePoleNumber(); p++) {
////        getDrudeResidue(p) *= Constants::eps0;
////    }
////    // Creates list of elements containing dispersive material.
////    vector<size_t> rpList;
////    rpList.reserve(cells.getLocalSize());
////    for (size_t e = 0; e < cells.getLocalSize(); e++) {
////        size_t id = cells.getIdOfRelPos(e);
////        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
////        if (cell->material->getId() == getId()) {
////            rpList.push_back(e);
////        }
////    }
////    nElem = rpList.size();
////    elem = new size_t[nElem];
////    for (size_t i = 0; i < nElem; i++) {
////        elem[i] = rpList[i];
////    }
////    // Usual poles.
////    dof = np * nElem * getPoleNumber();
////    P.setSize(dof);
////    P.setToZero();
////    rhsP.setSize(dof);
////    rhsP.setToZero();
////    resP.setSize(dof);
////    resP.setToZero();
////    // Drude poles.
////    drudeDof = np * nElem * getDrudePoleNumber();
////    J.setSize(drudeDof);
////    J.setToZero();
////    rhsJ.setSize(drudeDof);
////    rhsJ.setToZero();
////    resJ.setSize(drudeDof);
////    resJ.setToZero();
//}
//
