#include "Volume.h"

//DGDispersiveVolumic::DGDispersiveVolumic(const PMVolumeDispersive& mat) :
//                 PMVolumeDispersive(mat) 
// {
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
//void DGDispersiveVolumic::updateWithRes(
//        const size_t e1,
//        const size_t e2,
//        const Math::Real rkb) {
//    size_t i, e;
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            P.set(x)[i] += resP(x)[i] * rkb;
//            P.set(y)[i] += resP(y)[i] * rkb;
//            P.set(z)[i] += resP(z)[i] * rkb;
//        }
//    }
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
//        const Math::Real rka,
//        const Math::Real dt) {
//    size_t i, e;
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
//
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
//void DGDispersiveVolumic::computeRHSMagneticPolarizationCurrents(
//        const FieldR3& E, const size_t e1, const size_t e2) {
//}
//
//void DGDispersiveVolumic::addJumps(FieldR3& dE,
//        FieldR3& dH, FieldR3& E, FieldR3& H,
//        const size_t e1, const size_t e2) {
//}

