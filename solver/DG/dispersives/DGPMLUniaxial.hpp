// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
/*
 * SolverPMLUniaxial.cpp
 *
 *  Created on: Aug 2, 2013
 *      Author: luis
 */

#include "DGPMLUniaxial.h"

template<Int D>
DGPMLUniaxial<D>::DGPMLUniaxial(
        const PMVolumePML& mat,
        const CellGroup& cells,
        const bool useConductivity,
        const Real conductivity) :
DGPML(mat) {
    J.set(dof, 0.0);
    resJ.set(dof, 0.0);
    rhsJ.set(dof, 0.0);
    M.set(dof, 0.0);
    resM.set(dof, 0.0);
    rhsM.set(dof, 0.0);
    assert(check());
}

template<Int D>
DGPMLUniaxial<D>::~DGPMLUniaxial() {

}

template<Int D>
void DGPMLUniaxial<D>::addRHSToRes(
        const UInt e1, const UInt e2,
        const Real rka, const Real dt) {
//    UInt i,e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,e)
//#endif
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            resJ[i] *= rka;
//            resJ[i] += rhsJ[i] * dt;
//            resM[i] *= rka;
//            resM[i] += rhsM[i] * dt;
//        }
//    }
}

template<Int D>
void DGPMLUniaxial<D>::computeRHSElectric(
        FieldR3& rhsE, const FieldR3& E,
        const UInt e1, const UInt e2) const {
//    if (useConstantConductivity) {
//        UInt i, j, e, n;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n)
//#endif
//        for (i = 0; i < dof; i++) {
//            e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                n = i % np;
//                j = e * np + n ;
//                rhsE(dir1)[j] += (Constants::eps0*sig) * E(dir1)[j] - Constants::eps0 * J[i];
//                rhsE(dir2)[j] -= (Constants::eps0*sig) * E(dir2)[j];
//                rhsE(dir3)[j] -= (Constants::eps0*sig) * E(dir3)[j];
//            }
//        }
//    } else {
//        UInt i,j,e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e)
//#endif
//        for (e = 0; e < nElem; e++) {
//            if (e1 <= elem[e] && elem[e] < e2) {
//                i = e * np;
//                j = elem[e] * np;
//                //rhsE1[j] += (Constants::eps0*sig1) * E1[j] - Constants::eps0 * J[i];
//                add_am_v_prod<Real,np,np>(&rhsE.set(dir1)[j], sig1[e], &E(dir1)[j], Constants::eps0);
//                sub_a_v_prod<Real,np>(&rhsE(dir1)[j], &J[i], Constants::eps0);
//                //rhsE2[j] -= (Constants::eps0*sig1) * E2[j];
//                sub_am_v_prod<Real,np,np>(&rhsE(dir2)[j], sig1[e], &E(dir2)[j], Constants::eps0);
//                //rhsE3[j] -= (Constants::eps0*sig1) * E3[j];
//                sub_am_v_prod<Real,np,np>(&rhsE(dir3)[j], sig1[e], &E(dir3)[j], Constants::eps0);
//            }
//        }
//    }
}

template<Int D>
void DGPMLUniaxial<D>::computeRHSMagnetic(
        FieldR3& rhsH, const FieldR3& H,
        const UInt e1, const UInt e2) const {
//    if (useConstantConductivity) {
//        UInt i, j, e, n;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n)
//#endif
//        for (i = 0; i < dof; i++) {
//            e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                n = i % np;
//                j = e * np + n ;
//                rhsH(dir1)[j] += (Constants::mu0*sig) * H(dir1)[j] - Constants::mu0 * M[i];
//                rhsH(dir2)[j] -= (Constants::mu0*sig) * H(dir2)[j];
//                rhsH(dir3)[j] -= (Constants::mu0*sig) * H(dir3)[j];
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
//                //rhsH1[j] += (Constants::mu0*sigma1) * H1[j] - Constants::mu0 * M[i];
//                add_am_v_prod<Real,np,np>(&rhsH(dir1)[j], sig1[e], &H(dir1)[j], Constants::mu0);
//                sub_a_v_prod<Real,np>(&rhsH(dir1)[j], &M[i], Constants::mu0);
//                //rhsH2[j] -= (Constants::mu0*sigma1) * H2[j];
//                sub_am_v_prod<Real,np,np>(&rhsH(dir2)[j], sig1[e], &H(dir2)[j], Constants::mu0);
//                //rhsH3[j] -= (Constants::mu0*sigma1) * H3[j];
//                sub_am_v_prod<Real,np,np>(&rhsH(dir3)[j], sig1[e], &H(dir3)[j], Constants::mu0);
//            }
//        }
//    }
}

template<Int D>
void DGPMLUniaxial<D>::computeRHSElectricPolarizationCurrents(
        const FieldR3& E,
        const UInt e1, const UInt e2) {
//    if (useConstantConductivity) {
//        UInt i, j, e, n;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n)
//#endif
//        for (i = 0; i < dof; i++) {
//            e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                n = i % np;
//                j = e * np + n ;
//                rhsJ[i] = E(dir1)[j] * (sig*sig) - sig * J[i];
//            }
//        }
//    } else {
//        UInt i,j,e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e)
//#endif
//        for (e = 0; e < nElem; e++) {
//            if (e1 <= elem[e] && elem[e] < e2) {
//                i = e * np;
//                j = elem[e] * np;
//                //rhsJ[i] = E1[j] * (sig11) - sig1 * J[i];
//                m_v_prod<Real,np,np>(&rhsJ[i], sig11[e], &E(dir1)[j]);
//                sub_m_v_prod<Real,np,np>(&rhsJ[i], sig1[e], &J[i]);
//            }
//        }
//    }
}

template<Int D>
void DGPMLUniaxial<D>::computeRHSMagneticPolarizationCurrents(
        const FieldR3& H,
        const UInt e1, const UInt e2) {
//    if (useConstantConductivity) {
//        UInt i, j, e, n;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i,j,e,n)
//#endif
//        for (i = 0; i < dof; i++) {
//            e = elem[(i / np) % nElem];
//            if (e1 <= e && e < e2) {
//                n = i % np;
//                j = e * np + n ;
//                rhsM[i] = H(dir1)[j] * (sig*sig) - sig * M[i];
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
//                //rhsM[i] = H1[j] * (sigma1*sigma1) - sigma1 * M[i];
//                m_v_prod<Real,np,np>(&rhsM[i], sig11[e], &H(dir1)[j]);
//                sub_m_v_prod<Real,np,np>(&rhsM[i], sig1[e], &M[i]);
//            }
//        }
//    }
}

template<Int D>
bool DGPMLUniaxial<D>::check() const {
    bool sigInitialized = true;
    if (!useConstantConductivity) {
        sigInitialized &= (sig1 != NULL);
        sigInitialized &= (sig2 == NULL);
        sigInitialized &= (sig3 == NULL);
        sigInitialized &= (sig11 != NULL);
        sigInitialized &= (sig22 == NULL);
        sigInitialized &= (sig33 == NULL);
        sigInitialized &= (sig12 == NULL);
        sigInitialized &= (sig23 == NULL);
        sigInitialized &= (sig31 == NULL);
    }
    return sigInitialized;
}

template<Int D>
void DGPMLUniaxial<D>::updateWithRes(
        const UInt e1,
        const UInt e2,
        const Real rkb) {
    UInt i, e;
#ifdef SOLVER_USE_OPENMP
#pragma omp parallel for private(i, e)
#endif
    for (i = 0; i < dof; i++) {
        e = elem[(i / np) % nElem];
        if (e1 <= e && e < e2) {
            J[i] += resJ[i] * rkb;
            M[i] += resM[i] * rkb;
        }
    }
}
