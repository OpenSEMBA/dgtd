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
// * DGPMLTriaxial.cpp
// *
// *  Created on: Jun 22, 2015
// *      Author: luis
// */
//
//#include "DGPMLTriaxial.h"
//
//DGPMLTriaxial::DGPMLTriaxial(
//        const PMVolumePML& mat,
//        const CellGroup& cells,
//        const bool useConductivity,
//        const Real conductivity) :
//        DGPMLMultiaxial(mat, cells, useConductivity, conductivity) {
//}
//
//DGPMLTriaxial::~DGPMLTriaxial() {
//    // TODO Auto-generated destructor stub
//}
//
//void DGPMLTriaxial::computeRHSElectric(
//        FieldR3& rhsE,
//        const FieldR3& E,
//        const UInt e1, const UInt e2) const {
////    if (useConstantConductivity) {
////        UInt i, j, e, n;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i,j,e,n)
////#endif
////        for (i = 0; i < dof; i++) {
////            e = elem[(i / np) % nElem];
////            if (e1 <= e && e < e2) {
////                n = i % np;
////                j = e * np + n;
////                rhsE.set(0)[j] += - E(0)[j]*(Constants::eps0*sig) - J1[i]*Constants::eps0;
////                rhsE.set(0)[j] += - E(0)[j]*(Constants::eps0*sig) - J2[i]*Constants::eps0;
////                rhsE.set(0)[j] += - E(0)[j]*(Constants::eps0*sig) - J3[i]*Constants::eps0;
////            }
////        }
////    } else {
////        UInt i, j, e;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i,j,e)
////#endif
////        for (e = 0; e < nElem; e++) {
////            if (e1 <= elem[e] && e < elem[e]) {
////                i = e * np;
////                j = elem[e] * np;
////                //rhsEx[j] += - Ex[j]*Constants::eps0*(sig3+sig2-sig1) - J1[i]*Constants::eps0;
////                sub_am_v_prod<Real,np,np>(&rhsE.set(0)[j], sig3[e], &E(0)[j], Constants::eps0);
////                sub_am_v_prod<Real,np,np>(&rhsE.set(0)[j], sig2[e], &E(0)[j], Constants::eps0);
////                add_am_v_prod<Real,np,np>(&rhsE.set(0)[j], sig1[e], &E(0)[j], Constants::eps0);
////                sub_a_v_prod<Real,np>(&rhsE.set(0)[j], &J1[i], Constants::eps0);
////                //rhsEy[j] += - Ey[j]*Constants::eps0*(sig1+sig3-sig2) - J2[i]*Constants::eps0;
////                sub_am_v_prod<Real,np,np>(&rhsE.set(1)[j], sig1[e], &E(1)[j], Constants::eps0);
////                sub_am_v_prod<Real,np,np>(&rhsE.set(1)[j], sig3[e], &E(1)[j], Constants::eps0);
////                add_am_v_prod<Real,np,np>(&rhsE.set(1)[j], sig2[e], &E(1)[j], Constants::eps0);
////                sub_a_v_prod<Real,np>(&rhsE.set(1)[j], &J2[i], Constants::eps0);
////                //rhsEz[j] += - Ez[j]*Constants::eps0*(sig2+sig1-sig3) - J3[i]*Constants::eps0;
////                sub_am_v_prod<Real,np,np>(&rhsE.set(2)[j], sig2[e], &E(2)[j], Constants::eps0);
////                sub_am_v_prod<Real,np,np>(&rhsE.set(2)[j], sig1[e], &E(2)[j], Constants::eps0);
////                add_am_v_prod<Real,np,np>(&rhsE.set(2)[j], sig3[e], &E(2)[j], Constants::eps0);
////                sub_a_v_prod<Real,np>(&rhsE.set(2)[j], &J3[i], Constants::eps0);
////            }
////        }
////    }
//}
//
//void DGPMLTriaxial::computeRHSMagnetic(
//        FieldR3& rhsH,
//        const FieldR3& H, const UInt e1, const UInt e2) const {
////    if (useConstantConductivity) {
////        UInt i, j, e, n;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i,j,e,n)
////#endif
////        for (i = 0; i < dof; i++) {
////            e = elem[(i / np) % nElem];
////            if (e1 <= e && e < e2) {
////                n = i % np;
////                j = e * np + n;
////                rhsH.set(0)[j] += - H(0)[j]*(Constants::mu0*sig) - M1[i]*Constants::mu0;
////                rhsH.set(0)[j] += - H(0)[j]*(Constants::mu0*sig) - M2[i]*Constants::mu0;
////                rhsH.set(0)[j] += - H(0)[j]*(Constants::mu0*sig) - M3[i]*Constants::mu0;
////            }
////        }
////    } else {
////        UInt i, j, e;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i,j,e)
////#endif
////        for (e = 0; e < nElem; e++) {
////            if (e1 <= elem[e] && elem[e] < e2) {
////                i = e * np;
////                j = elem[e] * np;
////                //rhsHx[j] += - Hx[j]*Constants::mu0*(sig3+sig2-sig1) - M1[i]*Constants::mu0;
////                sub_am_v_prod<Real,np,np>(&rhsH.set(0)[j], sig3[e], &H(0)[j], Constants::mu0);
////                sub_am_v_prod<Real,np,np>(&rhsH.set(0)[j], sig2[e], &H(0)[j], Constants::mu0);
////                add_am_v_prod<Real,np,np>(&rhsH.set(0)[j], sig1[e], &H(0)[j], Constants::mu0);
////                sub_a_v_prod<Real,np>(&rhsH.set(0)[j], &M1[i], Constants::mu0);
////                //rhsHy[j] += - Hy[j]*Constants::mu0*(sig1+sig3-sig2) - M2[i]*Constants::mu0;
////                sub_am_v_prod<Real,np,np>(&rhsH.set(1)[j], sig1[e], &H(1)[j], Constants::mu0);
////                sub_am_v_prod<Real,np,np>(&rhsH.set(1)[j], sig3[e], &H(1)[j], Constants::mu0);
////                add_am_v_prod<Real,np,np>(&rhsH.set(1)[j], sig2[e], &H(1)[j], Constants::mu0);
////                sub_a_v_prod<Real,np>(&rhsH.set(1)[j], &M2[i], Constants::mu0);
////                //rhsHz[j] += - Hz[j]*Constants::mu0*(sig2+sig1-sig3) - M3[i]*Constants::mu0;
////                sub_am_v_prod<Real,np,np>(&rhsH.set(2)[j], sig2[e], &H(2)[j], Constants::mu0);
////                sub_am_v_prod<Real,np,np>(&rhsH.set(2)[j], sig1[e], &H(2)[j], Constants::mu0);
////                add_am_v_prod<Real,np,np>(&rhsH.set(2)[j], sig3[e], &H(2)[j], Constants::mu0);
////                sub_a_v_prod<Real,np>(&rhsH.set(2)[j], &M3[i], Constants::mu0);
////            }
////        }
////    }
//}
//
//void DGPMLTriaxial::computeRHSElectricPolarizationCurrents(
//        const FieldR3& E, const UInt e1, const UInt e2) {
////    if (useConstantConductivity) {
////        UInt i;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i)
////#endif
////        for (i = 0; i < dof; i++) {
////            UInt e = elem[(i / np) % nElem];
////            if (e1 <= e && e < e2) {
////                rhsJ1[i] = - J1[i] * sig;
////                rhsJ2[i] = - J2[i] * sig;
////                rhsJ3[i] = - J3[i] * sig;
////            }
////        }
////    } else {
////        UInt i, j, e;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i,j,e)
////#endif
////        for (e = 0; e < nElem; e++) {
////            if (e1 <= elem[e] && e < elem[e]) {
////                i = e * np;
////                j = elem[e] * np;
////                //rhsJ1[i] = sig23*Ex[j]-sig31*Ex[j]-sig12*Ex[j]+ sig11*Ex[j] - sig1*J1[i];
////                m_v_prod<Real,np,np>(&rhsJ1[i], sig23[e], &E(0)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ1[i], sig31[e], &E(0)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ1[i], sig12[e], &E(0)[j]);
////                add_m_v_prod<Real,np,np>(&rhsJ1[i], sig11[e], &E(0)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ1[i], sig1[e], &J1[i]);
////                //rhsJ2[i] = (sig1-sig2)*(sig3-sig2)*Ey[j] - sig2*J2[i];
////                m_v_prod<Real,np,np>(&rhsJ2[i], sig31[e], &E(1)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ2[i], sig12[e], &E(1)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ2[i], sig23[e], &E(1)[j]);
////                add_m_v_prod<Real,np,np>(&rhsJ2[i], sig22[e], &E(1)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ2[i], sig2[e], &J2[i]);
////                //rhsJ3[i] = (sig2-sig3)*(sig1-sig3)*E(2)[j] - sig3*J3[i];
////                m_v_prod<Real,np,np>(&rhsJ3[i], sig12[e], &E(2)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ3[i], sig23[e], &E(2)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ3[i], sig31[e], &E(2)[j]);
////                add_m_v_prod<Real,np,np>(&rhsJ3[i], sig33[e], &E(2)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsJ3[i], sig3[e], &J3[i]);
////            }
////        }
////    }
//}
//
//void DGPMLTriaxial::computeRHSMagneticPolarizationCurrents(
//        const FieldR3& H,
//        const UInt e1, const UInt e2) {
////    if (useConstantConductivity) {
////        UInt i;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i)
////#endif
////        for (i = 0; i < dof; i++) {
////            UInt e = elem[(i / np) % nElem];
////            if (e1 <= e && e < e2) {
////                rhsM1[i] = - M1[i]*sig;
////                rhsM2[i] = - M2[i]*sig;
////                rhsM3[i] = - M3[i]*sig;
////            }
////        }
////    } else {
////        UInt i, j, e;
////#ifdef SOLVER_USE_OPENMP
////#pragma omp parallel for private(i,j,e)
////#endif
////        for (e = 0; e < nElem; e++) {
////            if (e1 <= elem[e] && elem[e] < e2) {
////                i = e * np;
////                j = elem[e] * np;
////                //rhsM1[i] = (sig3 - sig1)*(sig2 - sig1)*Hx[j] - sig1*M1[i];
////                m_v_prod<Real,np,np>(&rhsM1[i], sig23[e], &H(0)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM1[i], sig31[e], &H(0)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM1[i], sig12[e], &H(0)[j]);
////                add_m_v_prod<Real,np,np>(&rhsM1[i], sig11[e], &H(0)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM1[i], sig1[e], &M1[i]);
////                //rhsM2[i] = (sig1 - sig2)*(sig3 - sig2)*Hy[j] - sig2*M2[i];
////                m_v_prod<Real,np,np>(&rhsM2[i], sig31[e], &H(1)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM2[i], sig12[e], &H(1)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM2[i], sig23[e], &H(1)[j]);
////                add_m_v_prod<Real,np,np>(&rhsM2[i], sig22[e], &H(1)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM2[i], sig2[e], &M2[i]);
////                //rhsM3[i] = (sig2 - sig3)*(sig1 - sig3)*Hz[j] - sig3*M3[i];
////                m_v_prod<Real,np,np>(&rhsM3[i], sig12[e], &H(2)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM3[i], sig23[e], &H(2)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM3[i], sig31[e], &H(2)[j]);
////                add_m_v_prod<Real,np,np>(&rhsM3[i], sig33[e], &H(2)[j]);
////                sub_m_v_prod<Real,np,np>(&rhsM3[i], sig3[e], &M3[i]);
////            }
////        }
////    }
//}
