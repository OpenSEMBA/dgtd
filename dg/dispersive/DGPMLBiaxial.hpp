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
// * DGPMLBiaxial.cpp
// *
// *  Created on: Jun 22, 2015
// *      Author: luis
// */
//
//#include "DGPMLBiaxial.h"
//
//DGPMLBiaxial::DGPMLBiaxial() {
//    // TODO Auto-generated constructor stub
//
//}
//
//DGPMLBiaxial::~DGPMLBiaxial() {
//    // TODO Auto-generated destructor stub
//}
//
//
//void
//DGPMLMultiaxial::internalBiaxialRHSElectricPolarizationCurrent(
//        const Real* E1, const Real* E2, const Real* E3, const UInt e1,
//        const UInt e2) {
//    if (useConstantConductivity) {
//        UInt i;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i)
//#endif
//        for (i = 0; i < dof; i++) {
//            UInt e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                rhsJ1[i] = - J1[i]*sig;
//                rhsJ2[i] = - J2[i]*sig;
//                rhsJ3[i] = J3[i] * (sig*sig);
//            }
//        }
//    } else {
//        UInt i, j, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e)
//#endif
//        for (e = 0; e < nElem; e++) {
//            if (e1 <= elem[e] && elem[e] < e2) {
//                i = e * np;
//                j = elem[e] * np;
//                //rhsJ1[i] = + sig11*E1[j] - sig12*E1[j] - sig1*J1[i];
//                m_v_prod<Real,np,np>(&rhsJ1[i], sig11[e], &E1[j]);
//                sub_m_v_prod<Real,np,np>(&rhsJ1[i], sig12[e], &E1[j]);
//                sub_m_v_prod<Real,np,np>(&rhsJ1[i], sig1[e], &J1[i]);
//                //rhsJ2[i] = - sig21*E2[j] + sig22*E2[j] - sig2*J2[i];
//                m_v_prod<Real,np,np>(&rhsJ2[i], sig22[e], &E2[j]);
//                sub_m_v_prod<Real,np,np>(&rhsJ2[i], sig12[e], &E2[j]);
//                sub_m_v_prod<Real,np,np>(&rhsJ2[i], sig2[e], &J2[i]);
//                //rhsJ3[i] = sig12*J3[i];
//                m_v_prod<Real,np,np>(&rhsJ3[i], sig12[e], &J3[i]);
//            }
//        }
//    }
//}
//
//void
//DGPMLMultiaxial::internalBiaxialRHSMagneticPolarizationCurrent(
//        const Real* H1, const Real* H2, const Real* H3,
//        const UInt e1, const UInt e2) {
//    if (useConstantConductivity) {
//        UInt i;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i)
//#endif
//        for (i = 0; i < dof; i++) {
//            UInt e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                rhsM1[i] = - M1[i]*sig;
//                rhsM2[i] = - M2[i]*sig;
//                rhsM3[i] = M3[i] * (sig*sig);
//            }
//        }
//    } else {
//        UInt i, j, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e)
//#endif
//        for (e = 0; e < nElem; e++) {
//            if (e1 <= e && e < e2) {
//                i = e * np;
//                j = elem[e] * np;
//                //rhsM1[i] = - sig12*H1[j] + sig11*H1[j] - sig1*M1[i];
//                m_v_prod<Real,np,np>(&rhsM1[i], sig11[e], &H1[j]);
//                sub_m_v_prod<Real,np,np>(&rhsM1[i], sig12[e], &H1[j]);
//                sub_m_v_prod<Real,np,np>(&rhsM1[i], sig1[e], &M1[i]);
//                //rhsM2[i] = - sig21*H2[j] + sig22*H2[j] - sig2*M2[i];
//                m_v_prod<Real,np,np>(&rhsM2[i], sig22[e], &H2[j]);
//                sub_m_v_prod<Real,np,np>(&rhsM2[i], sig12[e], &H2[j]);
//                sub_m_v_prod<Real,np,np>(&rhsM2[i], sig2[e], &M2[i]);
//                //rhsM3[i] = sig12*M3[i];
//                m_v_prod<Real,np,np>(&rhsM3[i], sig12[e], &M3[i]);
//            }
//        }
//    }
//}
//
//void
//DGPMLMultiaxial::internalBiaxialRHSMagnetic(
//        Real* rhsH1, Real* rhsH2, Real* rhsH3, const Real* H1, const Real* H2,
//        const Real* H3, const UInt e1, const UInt e2) const {
//    if (useConstantConductivity) {
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n)
//#endif
//        UInt i, j, e, n;
//        for (i = 0; i < dof; i++) {
//            e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                n = i % np;
//                j = e * np + n;
//                rhsH1[j] -= M1[i] * Constants::mu0;
//                rhsH2[j] -= M2[i] * Constants::mu0;
//                rhsH3[j] -= M3[i] * Constants::mu0;
//            }
//        }
//    } else {
//        UInt i, j, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e)
//#endif
//        for (e = 0; e < nElem; e++) {
//            if (e1 <= e && e < e2) {
//                i = e * np;
//                j = elem[e] * np;
//                //rhsM3[i] = sig12*M3[i];
//                m_v_prod<Real,np,np>(&rhsM3[i], sig12[e], &M3[i]);
//                //rhsH1[j] += - H1[j]*Constants::mu0*(sigma2-sigma1) - M1[i]*Constants::mu0;
//                sub_am_v_prod<Real,np,np>(&rhsH1[j], sig2[e], &H1[j], Constants::mu0);
//                add_am_v_prod<Real,np,np>(&rhsH1[j], sig1[e], &H1[j], Constants::mu0);
//                sub_a_v_prod<Real,np>(&rhsH1[j], &M1[i], Constants::mu0);
//                //rhsH2[j] += - H2[j]*Constants::mu0*(sigma1-sigma2) - M2[i]*Constants::mu0;
//                sub_am_v_prod<Real,np,np>(&rhsH2[j], sig1[e], &H2[j], Constants::mu0);
//                add_am_v_prod<Real,np,np>(&rhsH2[j], sig2[e], &H2[j], Constants::mu0);
//                sub_a_v_prod<Real,np>(&rhsH2[j], &M2[i], Constants::mu0);
//                //rhsH3[j] += - H3[j]*Constants::mu0*(sigma1+sigma2) - M3[i]*Constants::mu0;
//                sub_am_v_prod<Real,np,np>(&rhsH3[j], sig1[e], &H3[j], Constants::mu0);
//                sub_am_v_prod<Real,np,np>(&rhsH3[j], sig2[e], &H3[j], Constants::mu0);
//                sub_a_v_prod<Real,np>(&rhsH3[j], &J3[i], Constants::mu0);
//            }
//        }
//    }
//}
//
//
//void
//DGPMLMultiaxial::internalBiaxialRHSElectric(Real* rhsE1,
//        Real* rhsE2, Real* rhsE3, const Real* E1, const Real* E2,
//        const Real* E3, const UInt e1, const UInt e2) const {
//    if (useConstantConductivity) {
//        UInt i, j, e, n;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n)
//#endif
//        for (i = 0; i < dof; i++) {
//            e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                n = i % np;
//                j = e * np + n;
//                rhsE1[j] += - J1[i]*Constants::eps0;
//                rhsE2[j] += - J2[i]*Constants::eps0;
//                rhsE3[j] += - J3[i]*Constants::eps0;
//            }
//        }
//    } else {
//        UInt i, j, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e)
//#endif
//        for (e = 0; e < nElem; e++) {
//            if (e1 <= elem[e] && elem[e] < e2) {
//                i = e * np;
//                j = elem[e] * np;
//                //rhsE1[j] += - sig2*E1[j]*Constants::eps0 + sig1*E1[j]*Constants::eps0 - J1[i]*Constants::eps0;
//                sub_am_v_prod<Real,np,np>(&rhsE1[j], sig2[e], &E1[j], Constants::eps0);
//                add_am_v_prod<Real,np,np>(&rhsE1[j], sig1[e], &E1[j], Constants::eps0);
//                sub_a_v_prod<Real,np>(&rhsE1[j], &J1[i], Constants::eps0);
//                //rhsE2[j] += - E2[j]*Constants::eps0*(sig1-sig2) - J2[i]*Constants::eps0;
//                sub_am_v_prod<Real,np,np>(&rhsE2[j], sig1[e], &E2[j], Constants::eps0);
//                add_am_v_prod<Real,np,np>(&rhsE2[j], sig2[e], &E2[j], Constants::eps0);
//                sub_a_v_prod<Real,np>(&rhsE2[j], &J2[i], Constants::eps0);
//                //rhsE3[j] += - E3[j]*Constants::eps0*(sig1+sig2) - J3[i]*Constants::eps0;
//                sub_am_v_prod<Real,np,np>(&rhsE3[j], sig1[e], &E3[j], Constants::eps0);
//                sub_am_v_prod<Real,np,np>(&rhsE3[j], sig2[e], &E3[j], Constants::eps0);
//                sub_a_v_prod<Real,np>(&rhsE3[j], &J3[i], Constants::eps0);
//            }
//        }
//    }
//}
//
//DGPMLxy::DGPMLxy() {
//
//}
//
//DGPMLxy::DGPMLxy(
//        const  PMVolumePML& mat_,
//        const CellGroup& cells,
//        const bool useConductivity,
//        const Real conductivity) {
//    useConstantConductivity = useConductivity;
//    if (conductivity != 0.0) {
//        sig = conductivity;
//    }
//    initMultiaxial(mat_, cells);
//}
//
//DGPMLxy::~DGPMLxy() {
//
//}
//
//void
//DGPMLxy::computeRHSElectric(
//        FieldR3& rhs, const FieldR3& f,
//        const UInt e1, const UInt e2) const {
//    internalBiaxialRHSElectric(
//            rhs.set(x),rhs.set(y),rhs.set(z), f(x),f(y),f(z), e1,e2);
//}
//
//void
//DGPMLxy::computeRHSMagnetic(
//        FieldR3& rhs, const FieldR3& f,
//        const UInt e1, const UInt e2) const {
//    internalBiaxialRHSMagnetic(
//            rhs.set(x),rhs.set(y),rhs.set(z), f(x),f(y),f(z), e1,e2);
//}
//
//void
//DGPMLxy::computeRHSElectricPolarizationCurrents(
//        const FieldR3& f, const UInt e1, const UInt e2) {
//    internalBiaxialRHSElectricPolarizationCurrent(f(x),f(y),f(z), e1,e2);
//}
//
//void
//DGPMLxy::computeRHSMagneticPolarizationCurrents(
//        const FieldR3& f, const UInt e1, const UInt e2) {
//    internalBiaxialRHSMagneticPolarizationCurrent(f(x),f(y),f(z), e1,e2);
//}
