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
// * SolverSIBC.cpp
// *
// *  Created on: Jul 1, 2013
// *      Author: luis
// */
//
//#include "DGSIBC.h"
//
//DGSIBC::DGSIBC(
//        const PMSurfaceSIBC& mat_,
//        const CellGroup& cells,
//        Int*** map_, const Int vmapM_[faces][nfp],
//        Real*** ExP_, Real*** EyP_, Real*** EzP_,
//        Real*** HxP_, Real*** HyP_, Real*** HzP_) : PMSurfaceSIBC(mat_) {
////    map = map_;
////    ExP = ExP_;
////    EyP = EyP_;
////    EzP = EzP_;
////    HxP = HxP_;
////    HyP = HyP_;
////    HzP = HzP_;
////    for (size_t i = 0; i < faces; i++) {
////        for (size_t j = 0; j < nfp; j++) {
////            vmapM[i][j] = vmapM_[i][j];
////        }
////    }
////    // Inits sibc positions.
////    vector<pair<size_t, size_t> > efList0, efListD;
////    efList0.reserve(cells.getLocalSize());
////    efListD.reserve(cells.getLocalSize());
////    for (size_t i = 0; i < bc.size(); i++) {
////        size_t id = bc[i]->getCell()->getId();
////        size_t f = bc[i]->getFace();
////        if (!bc[i]->isBack()) {
////            pair<size_t,size_t> aux0(cells.getRelPosOfId(id), f);
////            efList0.push_back(aux0);
////        } else {
////            pair<size_t,size_t> auxD(cells.getRelPosOfId(id), f);
////            efListD.push_back(auxD);
////        }
////        // TODO Only works with flat faces.
////        assert(!bc[i]->getCell()->isCurvedFace(f));
////        // TODO Can't be a computational domain boundary.
////    }
////    nE0 = efList0.size();
////    elem0 = new size_t[nE0];
////    face0 = new size_t[nE0];
////    for (size_t i = 0; i < nE0; i++) {
////        elem0[i] = efList0[i].first;
////        face0[i] = efList0[i].second;
////    }
////    nED = efListD.size();
////    elemD = new size_t[nED];
////    faceD = new size_t[nED];
////    for (size_t i = 0; i < nED; i++) {
////        elemD[i] = efListD[i].first;
////        faceD[i] = efListD[i].second;
////    }
////    // Normal vectors.
////    n0 = new CVecR3[nE0];
////    for (size_t k = 0; k < nE0; k++) {
////        n0[k] = bc[k]->getCell()->sideNormal(face0[k]);
////    }
////    nD = new CVecR3[nED];
////    for (size_t k = 0; k < nED; k++) {
////        nD[k] = bc[k]->getCell()->sideNormal(faceD[k]);
////    }
////    // Allocates magnetic currents.
////    nP = getNumberOfPoles();
////    Q0 = new CVecR3*[nP];
////    rhsQ0 = new CVecR3*[nP];
////    resQ0 = new CVecR3*[nP];
////    QD = new CVecR3*[nP];
////    rhsQD = new CVecR3*[nP];
////    resQD = new CVecR3*[nP];
////    for (size_t p = 0; p < nP; p++) {
////        Q0[p] = new CVecR3[nE0*nfp];
////        rhsQ0[p] = new CVecR3[nE0*nfp];
////        resQ0[p] = new CVecR3[nE0*nfp];
////    }
////    for (size_t p = 0; p < nP; p++) {
////        QD[p] = new CVecR3[nED*nfp];
////        rhsQD[p] = new CVecR3[nED*nfp];
////        resQD[p] = new CVecR3[nED*nfp];
////    }
////    // Polarization fields.
////    E0 = new CVecR3[nE0*nfp];
////    ED = new CVecR3[nED*nfp];
//}
//
//
//DGSIBC::~DGSIBC() {
//
//}
//
//void
//DGSIBC::computeRHSMagneticPolarizationCurrents(
//        const FieldR3& H,
//        const size_t e1, const size_t e2) {
//    size_t e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(e)
//#endif
//    for (e = 0; e < nE0; e++) {
//        size_t e0 = elem0[e]; // Element-face pair.
//        size_t f0 = face0[e];
//        size_t eD = elemD[e];
//        size_t fD = faceD[e];
//        if (e1 <= e0 && e0 < e2) {
//            for (size_t p = 0; p < nP; p++) {
//                size_t i = e * nfp;
//                for (size_t j = 0; j < nfp; j++) {
//                    Int vM = e0 * np + vmapM[f0][j];
//                    const CVecR3 H0(H(x)[vM],H(y)[vM],H(z)[vM]);
//                    size_t vP = map[e0][f0][j]; // Field coeff pos.
//                    const CVecR3 HD(HxP[eD][fD][vP], HyP[eD][fD][vP], HzP[eD][fD][vP]);
//                    rhsQ0[p][i] =  (Q0[p][i] * pole_[p])
//                                           + (n0[e] ^ H0) * Z_[p](0,0) + (n0[e] ^ HD) * Z_[p](0,1);
//                    i++;
//                }
//            }
//        }
//    }
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(e)
//#endif
//    for (e = 0; e < nED; e++) {
//        size_t e0 = elem0[e]; // Element-face pair.
//        size_t f0 = face0[e];
//        size_t eD = elemD[e];
//        size_t fD = faceD[e];
//        if (e1 <= eD && eD < e2) {
//            for (size_t p = 0; p < nP; p++) {
//                size_t i = e * nfp; // Node number.
//                for (size_t j = 0; j < nfp; j++) {
//                    size_t vM = eD * np + vmapM[fD][j]; // Field coeff pos.
//                    const CVecR3
//                    HD(H(x)[vM],H(y)[vM],H(z)[vM]);
//                    size_t vP = map[eD][fD][j]; // Field coeff pos.
//                    const CVecR3
//                    H0(HxP[e0][f0][vP], HyP[e0][f0][vP], HzP[e0][f0][vP]);
//                    rhsQD[p][i] = QD[p][i] * pole_[p]
//                                                   - (nD[e] ^ H0) * Z_[p](1,0) - (nD[e] ^ HD) * Z_[p](1,1);
//                    i++;
//                }
//            }
//        }
//    }
//}
//
//void
//DGSIBC::computeRHSMagnetic(
//        FieldR3& rhsH,
//        const FieldR3& H,
//        const size_t e1, const size_t e2) const {
//
//}
//
//void
//DGSIBC::computePolarizationFields(
//        const Real* Hx, const Real* Hy,	const Real* Hz,
//        const size_t e1, const size_t e2) {
//    for (size_t e = 0; e < nE0; e++) {
//        size_t e0 = elem0[e];
//        size_t f0 = face0[e];
//        size_t eD = elemD[e];
//        size_t fD = faceD[e];
//        if (e1 <= e0 && e0 < e2) {
//            for (size_t j = 0; j < nfp; j++) {
//                Int vM = e0 * np + vmapM[f0][j];
//                CVecR3 H0(Hx[vM],Hy[vM],Hz[vM]);
//                size_t vP = map[eD][fD][j]; // Field coeff pos.
//                CVecR3
//                HD(HxP[eD][fD][vP], HyP[eD][fD][vP], HzP[eD][fD][vP]);
//                size_t i = e * nfp;
//#				warning "Ignoring Zinfinite"
//                //				E0[i] = (n0[e] ^ H0) * mat.Zinfinite[0]
//                //				 + (n0[e] ^ HD) * mat.Zinfinite[1];
//                for (size_t p = 0; p < nP; p++) {
//                    E0[i] += Q0[p][i];
//                }
//                i++;
//            }
//        }
//    }
//    for (size_t e = 0; e < nED; e++) {
//        size_t e0 = elem0[e];
//        size_t f0 = face0[e];
//        size_t eD = elemD[e];
//        size_t fD = faceD[e];
//        if (e1 <= e0 && e0 < e2) {
//            for (size_t j = 0; j < nfp; j++) {
//                Int vM = eD * np + vmapM[fD][j];
//                CVecR3 H0(Hx[vM],Hy[vM],Hz[vM]);
//                size_t vP = map[e0][f0][j]; // Field coeff pos.
//                CVecR3
//                HD(HxP[e0][f0][vP], HyP[e0][f0][vP], HzP[e0][f0][vP]);
//                size_t i = e * nfp;
//#				warning "Ignoring Zinfinite"
//                //				ED[i] = - (nD[e] ^ H0) * mat.Zinfinite[2]
//                //				 - (nD[e] ^ HD) * mat.Zinfinite[3];
//                for (size_t p = 0; p < nP; p++) {
//                    ED[i] += QD[p][i];
//                }
//                i++;
//            }
//        }
//    }
//}
//
//void
//DGSIBC::addJumps(
//        FieldR3& dE, FieldR3& dH,
//        FieldR3& E, FieldR3& H,
//        const size_t e1, const size_t e2) {
//    // Updates E0 and ED.
//    computePolarizationFields(H(x),H(y),H(z),e1,e2);
//    // Updates jumps with polarized fields.
//    size_t i, f, e, j, n;
//    for (size_t lE = 0; lE < nE0; lE++) {
//        f = face0[lE];
//        e = elem0[lE];
//        if (e >= e1 && e < e2) {
//            n = (e * faces + f) * nfp;
//            i = lE * nfp;
//            for (j = 0; j < nfp; j++) {
//                Int vM = e * np + vmapM[f][j];
//                dE.set(x)[n] = 2.0 * (E(x)[vM] - E0[i](0));
//                dE.set(y)[n] = 2.0 * (E(y)[vM] - E0[i](1));
//                dE.set(z)[n] = 2.0 * (E(z)[vM] - E0[i](2));
//                dH.set(x)[n] = 0.0;
//                dH.set(y)[n] = 0.0;
//                dH.set(z)[n] = 0.0;
//                i++;
//                n++;
//            }
//        }
//    }
//    for (size_t lE = 0; lE < nED; lE++) {
//        f = faceD[lE];
//        e = elemD[lE];
//        if (e >= e1 && e < e2) {
//            n = (e * faces + f) * nfp;
//            i = nfp * lE;
//            for (j = 0; j < nfp; j++) {
//                Int vM = e * np + vmapM[f][j];
//                dE.set(x)[n] = 2.0 * (E(x)[vM] - ED[i](0));
//                dE.set(y)[n] = 2.0 * (E(y)[vM] - ED[i](1));
//                dE.set(z)[n] = 2.0 * (E(z)[vM] - ED[i](2));
//                dH.set(x)[n] = 0.0;
//                dH.set(y)[n] = 0.0;
//                dH.set(z)[n] = 0.0;
//                i++;
//                n++;
//            }
//        }
//    }
//}
//
//void
//DGSIBC::addRHSToRes(
//        const size_t e1, const size_t e2,
//        const Real rka,
//        const Real dt) {
//    for (size_t e = 0; e < nE0; e++) {
//        if (e >= e1 && e < e2) {
//            for (size_t p = 0; p < nP; p++) {
//                size_t i = e*nfp;
//                for (size_t j = 0; j < nfp; j++) {
//                    resQ0[p][i] = resQ0[p][i] * rka + rhsQ0[p][i] * dt;
//                    i++;
//                }
//            }
//        }
//    }
//    for (size_t e = 0; e < nED; e++) {
//        if (e >= e1 && e < e2) {
//            for (size_t p = 0; p < nP; p++) {
//                size_t i = e*nfp;
//                for (size_t j = 0; j < nfp; j++) {
//                    resQD[p][i] = resQD[p][i] * rka + rhsQD[p][i] * dt;
//                    i++;
//                }
//            }
//        }
//    }
//}
//
//void
//DGSIBC::updateWithRes(
//        const size_t e1, const size_t e2, const Real rkb) {
//    for (size_t e = 0; e < nE0; e++) {
//        if (e >= e1 && e < e2) {
//            for (size_t p = 0; p < nP; p++) {
//                size_t i = e*nfp;
//                for (size_t j = 0; j < nfp; j++) {
//                    Q0[p][i] += resQ0[p][i] * rkb;
//                }
//            }
//        }
//    }
//    for (size_t e = 0; e < nED; e++) {
//        if (e >= e1 && e < e2) {
//            for (size_t p = 0; p < nP; p++) {
//                size_t i = e*nfp;
//                for (size_t j = 0; j < nfp; j++) {
//                    QD[p][i] += resQD[p][i] * rkb;
//                }
//            }
//        }
//    }
//}
//
//void
//DGSIBC::computeRHSElectric(
//        FieldR3& rhsE,
//        const FieldR3& E,
//        const size_t e1, const size_t e2) const {
//}
//
//void
//DGSIBC::computeRHSElectricPolarizationCurrents(
//        const FieldR3& E,
//        const size_t e1, const size_t e2) {
//}
