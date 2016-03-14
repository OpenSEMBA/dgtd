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
 * SolverExplicit.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: luis
 */

#include "../dg/DGExplicit.h"

DGExplicit::DGExplicit(
        const Mesh::Volume& mesh,
        const PMGroup& pMGroup,
        const EMSourceGroup& emSources,
        const OptionsSolverDGTD& options,
        Comm* comm) {
    CellGroup cells(mesh, pMGroup);
    init(options, pMGroup, cells, comm);
    allocateRHSAndJumps();
    if (options.isUseLTS()) {
        allocateFieldsForLTS();
    }
    if (emSources.size() != 0) {
        E.setAll((Real) 0.0);
        H.setAll((Real) 0.0);
    } else {
        setFieldsToRandom();
        cout<< ">> No EM Excitations were detected <<" << endl;
        cout<< ">> A random field is being used  <<" << endl;
    }
    allocateMaps();
    deduplicateVMaps(cells);

    Connectivities map(mesh.elems());
    BCGroup bc(mesh, emSources, pMGroup, cells, map);
    assignPointersToNeighbours(cells, map, mesh);
    buildEMSources(emSources, bc, map, cells);
    BCToLocalArray(bc, cells, map);
    buildScalingFactors(cells, map);
}

DGExplicit::~DGExplicit() {

}

size_t DGExplicit::getFieldDOFs() {
    return (nK * np * 3);
}

const FieldR3& DGExplicit::getRHSElectric() const {
    return rhsE;
}

const FieldR3& DGExplicit::getRHSMagnetic() const {
    return rhsH;
}

void DGExplicit:: printInfo() const {
    cout<< " --- SolverExplicit info ---" << endl;
    cout<< "- Boundary conditions: " << endl;
    cout<< "# SMA: " << nSMA << endl;
    cout<< "# PEC: " << nPEC << endl;
    cout<< "# PMC: " << nPMC << endl;
    cout<< "# Of Curved Faces: " << nCurvedFaces << endl;
}

void DGExplicit::computePolarizationCurrentsRHS(
        const size_t e1,
        const size_t e2) {
    computePolarizationCurrentsRHSElectric(e1,e2);
    computePolarizationCurrentsRHSMagnetic(e1,e2);
}

void DGExplicit::computePolarizationCurrentsRHSElectric(
        const size_t e1, const size_t e2) {
    for (size_t d = 0; d < dispersive.size(); d++) {
        dispersive[d]->computeRHSElectricPolarizationCurrents(E, e1, e2);
        dispersive[d]->computeRHSElectric(rhsE, E, e1, e2);
    }
}

void DGExplicit::computePolarizationCurrentsRHSMagnetic(
        const size_t e1, const size_t e2) {
    for (size_t d = 0; d < dispersive.size(); d++) {
        dispersive[d]->computeRHSMagneticPolarizationCurrents(H, e1, e2);
        dispersive[d]->computeRHSMagnetic(rhsH, H, e1, e2);
    }
}

void DGExplicit::computeRHS(
        const size_t e1,
        const size_t e2,
        const Real localTime,
        const Real minDT) {
    computeRHSElectric(e1,e2, localTime,minDT);
    computeRHSMagnetic(e1,e2, localTime,minDT);
}

void DGExplicit::computeRHSElectric(
        const size_t e1,
        const size_t e2,
        const Real localTime,
        const Real minDT) {
    computeCurlsInRHSElectric(e1,e2);
    computeJumps(e1,e2, localTime,minDT);
    addFluxesToRHSElectric(e1,e2);
    size_t i, j, e;
#	pragma omp parallel for private(i,j,e)
    for (e = e1; e < e2; e++) {
        i = e * np;
        j = (e + 1) * np;
        rhsE.prod(i, j, oneOverEps[e]);
    }
}

void DGExplicit::computeRHSMagnetic(
        const size_t e1,
        const size_t e2,
        const Real localTime,
        const Real minDT) {
    computeCurlsInRHSMagnetic(e1,e2);
    computeJumps(e1,e2,localTime,minDT);
    addFluxesToRHSMagnetic(e1,e2);
    size_t i, j, e;
#pragma omp parallel for private(i,j,e)
    for (e = e1; e < e2; e++) {
        i = e * np;
        j = (e + 1) * np;
        rhsH.prod(i, j, oneOverMu[e]);
    }
}

void DGExplicit::computeCurlsInRHSElectric(
        const size_t e1,
        const size_t e2) {
    size_t i, e;
#pragma omp parallel for private(e,i)
    for (e = e1; e < e2; e++) {
        // i: Beginning of element field. [0, (nK-1)*np]
        i = e * np;
        // rhsE(x) = + Cy * H(z) - Cz * H(y)
        m_v_prod<Real,np,np>(&rhsE.set(x)[i], Cy[e], &H(z)[i]);
        sub_m_v_prod<Real,np,np>(&rhsE.set(x)[i], Cz[e], &H(y)[i]);
        // rhsE(y) = + Cz * H(x) - Cx * H(z)
        m_v_prod<Real,np,np>(&rhsE.set(y)[i], Cz[e], &H(x)[i]);
        sub_m_v_prod<Real,np,np>(&rhsE.set(y)[i], Cx[e], &H(z)[i]);
        // rhsE(z) = + Cx * H(y) - Cy * H(x)
        m_v_prod<Real,np,np>(&rhsE.set(z)[i], Cx[e], &H(y)[i]);
        sub_m_v_prod<Real,np,np>(&rhsE.set(z)[i], Cy[e], &H(x)[i]);
    }
}

void DGExplicit::computeCurlsInRHSMagnetic(
        const size_t e1,
        const size_t e2) {
    size_t i, e;
#pragma omp parallel for private(e,i)
    for (e = e1; e < e2; e++) {
        // i: Beginning of element field. [0, (nK-1)*np]
        i = e * np;
        // rhsH(x) = - Cy * E(z) + Cz * E(y)
        am_v_prod<Real,np,np>(&rhsH.set(x)[i], Cy[e], &E(z)[i], -1.0);
        add_m_v_prod<Real,np,np>(&rhsH.set(x)[i], Cz[e], &E(y)[i]);
        // rhsH(y) = - Cz * E(x) + Cx * E(z)
        am_v_prod<Real,np,np>(&rhsH.set(y)[i], Cz[e], &E(x)[i], -1.0);
        add_m_v_prod<Real,np,np>(&rhsH.set(y)[i], Cx[e], &E(z)[i]);
        // rhsH(z) = - Cx * E(y) + Cy * E(x)
        am_v_prod<Real,np,np>(&rhsH.set(z)[i], Cx[e], &E(y)[i], -1.0);
        add_m_v_prod<Real,np,np>(&rhsH.set(z)[i], Cy[e], &E(x)[i]);
    }
}
void DGExplicit::computeJumps(
        const size_t e1,
        const size_t e2,
        const Real localTime,
        const Real minDT) {
    if (comm->getNumberOfTasks() > 1) {
        comm->syncNeighbourFields(
                nE.set(x), nE.set(y), nE.set(z), nH.set(x), nH.set(y), nH.set(z),
                E(x),E(y),E(z),H(x),H(y),H(z));
    }
    size_t b, i, j, k, f, e;
    size_t vM, vP;
#pragma omp parallel for private(e,k,i,f,j,vM,vP)
    for (e = e1; e < e2; e++) {
        k = e * np;   // Beginning of element field. [0, (nK-1)*np]
        i = e * nfp * faces;
        for (size_t f = 0; f < faces; f++) {
            for (size_t j = 0; j < nfp; j++) {
                vM = k + vmapM[f][j]; // Local field pos.
                vP = map_[e][f][j]; // Neigh field pos.
                dE.set(x)[i] = E(x)[vM] - ExP[e][f][vP];
                dE.set(y)[i] = E(y)[vM] - EyP[e][f][vP];
                dE.set(z)[i] = E(z)[vM] - EzP[e][f][vP];
                dH.set(x)[i] = H(x)[vM] - HxP[e][f][vP];
                dH.set(y)[i] = H(y)[vM] - HyP[e][f][vP];
                dH.set(z)[i] = H(z)[vM] - HzP[e][f][vP];
                i++;
            }
        }
    }
    // SMA
#pragma omp parallel for private(b,e,k,i,f,j,vM,vP)
    for (b = 0; b < nSMA; b++) {
        f = SMAf[b];
        e = SMAe[b];
        if (e >= e1 && e < e2) {
            i = e * nfp * faces + f * nfp;
            k = e * np;
            for (j = 0; j < nfp; j++) {
                vM = vmapM[f][j];
                dE.set(x)[i + j] = E(x)[k + vM];
                dE.set(y)[i + j] = E(y)[k + vM];
                dE.set(z)[i + j] = E(z)[k + vM];
                dH.set(x)[i + j] = H(x)[k + vM];
                dH.set(y)[i + j] = H(y)[k + vM];
                dH.set(z)[i + j] = H(z)[k + vM];
            }
        }
    }
    // PEC
#pragma omp parallel for private(b,e,k,i,f,j,vM,vP)
    for (b = 0; b < nPEC; b++) {
        f = PECf[b];
        e = PECe[b];
        if (e >= e1 && e < e2) {
            i = e * nfp * faces + f * nfp;
            k = e * np;
            for (j = 0; j < nfp; j++) {
                vM = vmapM[f][j];
                dE.set(x)[i + j] = 2.0 * E(x)[k + vM];
                dE.set(y)[i + j] = 2.0 * E(y)[k + vM];
                dE.set(z)[i + j] = 2.0 * E(z)[k + vM];
                dH.set(x)[i + j] = 0.0;
                dH.set(y)[i + j] = 0.0;
                dH.set(z)[i + j] = 0.0;
            }
        }
    }
    // PMC
#pragma omp parallel for private(b,e,k,i,f,j,vM,vP)
    for (b = 0; b < nPMC; b++) {
        f = PMCf[b];
        e = PMCe[b];
        if (e >= e1 && e < e2) {
            i = e * nfp * faces + f * nfp;
            k = e * np;
            for (j = 0; j < nfp; j++) {
                vM = vmapM[f][j];
                dE.set(x)[i + j] = 0.0;
                dE.set(y)[i + j] = 0.0;
                dE.set(z)[i + j] = 0.0;
                dH.set(x)[i + j] = 2.0 * H(x)[k + vM];
                dH.set(y)[i + j] = 2.0 * H(y)[k + vM];
                dH.set(z)[i + j] = 2.0 * H(z)[k + vM];
            }
        }
    }
    // Computes E(x)citations.
    for (b = 0; b < source.size(); b++) {
        source[b]->computeExcitation(localTime,minDT);
        source[b]->addJumps(e1,e2);
    }
    // Computes contributions of polarization currents to jumps
    for (b = 0; b < dispersive.size(); b++) {
        dispersive[b]->addJumps(dE, dH, E, H, e1,e2);
    }
}

void DGExplicit::addFluxesToRHSElectric(
        const size_t e1,
        const size_t e2) {
    static const bool useResForUpw = false;
    addFluxesToRHSElectric(e1,e2,useResForUpw);
}

void DGExplicit::addFluxesToRHSMagnetic(
        const size_t e1,
        const size_t e2) {
    static const bool useResForUpw = false;
    addFluxesToRHSMagnetic(e1,e2,useResForUpw);
}

void DGExplicit::addFluxesToRHSElectric(
        const size_t e1,
        const size_t e2,
        const bool useResForUpw) {
    addStraightFluxesToRHSElectric(e1,e2,useResForUpw);
    addCurvedFluxesToRHSElectric(e1,e2,useResForUpw);
}

void DGExplicit::addFluxesToRHSMagnetic(
        const size_t e1,
        const size_t e2,
        const bool useResForUpw) {
    addStraightFluxesToRHSMagnetic(e1,e2,useResForUpw);
    addCurvedFluxesToRHSMagnetic(e1,e2,useResForUpw);
}

void DGExplicit::addStraightFluxesToRHSElectric(
        const size_t e1,
        const size_t e2,
        const bool useResForUpw) {
    size_t i,j, k, f, e;
    Real fx[nfpfaces], fy[nfpfaces], fz[nfpfaces];
    if (upwinding == 0.0) {
        // ---------- Centred flux ------------------------------------
#ifdef SOLVER_USE_OPENMP
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
#endif
        for (e = e1; e < e2; e++) {
            i = e * nfpfaces;
            for (k = 0; k < nfpfaces; k++) {
                f = k / nfp;
                j = e * faces + f;
                // - n ^ dH ---------------------------------------------------
                fx[k] = - nImp(y)[j] * dH(z)[i] + nImp(z)[j] * dH(y)[i];
                fy[k] = - nImp(z)[j] * dH(x)[i] + nImp(x)[j] * dH(z)[i];
                fz[k] = - nImp(x)[j] * dH(y)[i] + nImp(y)[j] * dH(x)[i];
                i++;
            }
            i = e * np;
            add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(x)[i], LIFT, fx);
            add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(y)[i], LIFT, fy);
            add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(z)[i], LIFT, fz);
        }
    } else if (upwinding == 1.0) {
        // ---------- Upwind flux -------------------------------------
        if (useResForUpw) {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = - nImp(y)[j]  * dH(z)[i] + nImp(z)[j]  * dH(y)[i]
                                                                           + cnImp(x)[j] * dresE(y)[i] + cnImp(z)[j] * dresE(z)[i]
                                                                                                                                - rnImp(x)[j] * dresE(x)[i];
                    fy[k] = - nImp(z)[j]  * dH(x)[i] + nImp(x)[j]  * dH(z)[i]
                                                                           + cnImp(y)[j] * dresE(z)[i] + cnImp(x)[j] * dresE(x)[i]
                                                                                                                                - rnImp(y)[j] * dresE(y)[i];
                    fz[k] = - nImp(x)[j]  * dH(y)[i] + nImp(y)[j]  * dH(x)[i]
                                                                           + cnImp(z)[j] * dresE(x)[i] + cnImp(y)[j] * dresE(y)[i]
                                                                                                                                - rnImp(z)[j] * dresE(z)[i];
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(z)[i], LIFT, fz);
            }
        } else {
            // >>>>>>> This one is the most frequently used <<<<<<<<<
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = - nImp(y)[j]  * dH(z)[i] + nImp(z)[j]  * dH(y)[i]
                                                                           + cnImp(x)[j] * dE(y)[i] + cnImp(z)[j] * dE(z)[i]
                                                                                                                          - rnImp(x)[j] * dE(x)[i];
                    fy[k] = - nImp(z)[j]  * dH(x)[i] + nImp(x)[j]  * dH(z)[i]
                                                                           + cnImp(y)[j] * dE(z)[i] + cnImp(x)[j] * dE(x)[i]
                                                                                                                          - rnImp(y)[j] * dE(y)[i];
                    fz[k] = - nImp(x)[j]  * dH(y)[i] + nImp(y)[j]  * dH(x)[i]
                                                                           + cnImp(z)[j] * dE(x)[i] + cnImp(y)[j] * dE(y)[i]
                                                                                                                          - rnImp(z)[j] * dE(z)[i];
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(z)[i], LIFT, fz);
            }
            // >>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<
        }
    } else {
        // ---------- Penalized flux ----------------------------------
        if (useResForUpw) {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = - nImp(y)[j]  * dH(z)[i] + nImp(z)[j]  * dH(y)[i]
                                                                           + (cnImp(x)[j] * dresE(y)[i] + cnImp(z)[j] * dresE(z)[i]
                                                                                                                                 - rnImp(x)[j] * dresE(x)[i]) * upwinding;
                    fy[k] = - nImp(z)[j]  * dH(x)[i] + nImp(x)[j]  * dH(z)[i]
                                                                           + (cnImp(y)[j] * dresE(z)[i] + cnImp(x)[j] * dresE(x)[i]
                                                                                                                                 - rnImp(y)[j] * dresE(y)[i]) * upwinding;
                    fz[k] = - nImp(x)[j]  * dH(y)[i] + nImp(y)[j]  * dH(x)[i]
                                                                           + (cnImp(z)[j] * dresE(x)[i] + cnImp(y)[j] * dresE(y)[i]
                                                                                                                                 - rnImp(z)[j] * dresE(z)[i]) * upwinding;
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(z)[i], LIFT, fz);
            }
        } else {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = - nImp(y)[j]  * dH(z)[i] + nImp(z)[j]  * dH(y)[i]
                                                                           + (cnImp(x)[j] * dE(y)[i] + cnImp(z)[j] * dE(z)[i]
                                                                                                                           - rnImp(x)[j] * dE(x)[i]) * upwinding;
                    fy[k] = - nImp(z)[j]  * dH(x)[i] + nImp(x)[j]  * dH(z)[i]
                                                                           + (cnImp(y)[j] * dE(z)[i] + cnImp(x)[j] * dE(x)[i]
                                                                                                                           - rnImp(y)[j] * dE(y)[i]) * upwinding;
                    fz[k] = - nImp(x)[j]  * dH(y)[i] + nImp(y)[j]  * dH(x)[i]
                                                                           + (cnImp(z)[j] * dE(x)[i] + cnImp(y)[j] * dE(y)[i]
                                                                                                                           - rnImp(z)[j] * dE(z)[i]) * upwinding;
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(x)[i],LIFT,fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(y)[i],LIFT,fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsE.set(z)[i],LIFT,fz);
            }
        }
    }
}

void DGExplicit::addStraightFluxesToRHSMagnetic (
        const size_t e1,
        const size_t e2,
        const bool useResForUpw) {
    size_t i, j, f, e, k;
    Real fx[nfpfaces], fy[nfpfaces], fz[nfpfaces];
    if (upwinding == 0.0) {
        // ---------- Centred flux --------------------------------------------
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
        for (e = e1; e < e2; e++) {
            i = e * nfpfaces;
            for (k = 0; k < nfpfaces; k++) {
                f = k / nfp;
                j = e * faces + f;
                fx[k] = + nAdm(y)[j] * dE(z)[i] - nAdm(z)[j]  * dE(y)[i];
                fy[k] = + nAdm(z)[j] * dE(x)[i] - nAdm(x)[j]  * dE(z)[i];
                fz[k] = + nAdm(x)[j] * dE(y)[i] - nAdm(y)[j]  * dE(x)[i];
                i++;
            }
            i = e * np;
            add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(x)[i], LIFT, fx);
            add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(y)[i], LIFT, fy);
            add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(z)[i], LIFT, fz);
        }
    } else if (upwinding == 1.0) {
        // ---------- Upwind flux ---------------------------------------------
        if (useResForUpw) {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = + nAdm(y)[j]  * dE(z)[i] - nAdm(z)[j]  * dE(y)[i]
                                                                           + cnAdm(x)[j] * dresH(y)[i] + cnAdm(z)[j] * dresH(z)[i]
                                                                                                                                - rnAdm(x)[j] * dresH(x)[i];
                    fy[k] = + nAdm(z)[j]  * dE(x)[i] - nAdm(x)[j]  * dE(z)[i]
                                                                           + cnAdm(y)[j] * dresH(z)[i] + cnAdm(x)[j] * dresH(x)[i]
                                                                                                                                - rnAdm(y)[j] * dresH(y)[i];
                    fz[k] = + nAdm(x)[j]  * dE(y)[i]  - nAdm(y)[j] * dE(x)[i]
                                                                           + cnAdm(z)[j] * dresH(x)[i]  + cnAdm(y)[j] * dresH(y)[i]
                                                                                                                                 - rnAdm(z)[j] * dresH(z)[i];
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(z)[i], LIFT, fz);
            }
        } else {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = + nAdm(y)[j]  * dE(z)[i] - nAdm(z)[j]  * dE(y)[i]
                                                                           + cnAdm(x)[j] * dH(y)[i] + cnAdm(z)[j] * dH(z)[i]
                                                                                                                          - rnAdm(x)[j] * dH(x)[i];
                    fy[k] = + nAdm(z)[j]  * dE(x)[i] - nAdm(x)[j]  * dE(z)[i]
                                                                           + cnAdm(y)[j] * dH(z)[i] + cnAdm(x)[j] * dH(x)[i]
                                                                                                                          - rnAdm(y)[j] * dH(y)[i];
                    fz[k] = + nAdm(x)[j]  * dE(y)[i]  - nAdm(y)[j]  * dE(x)[i]
                                                                            + cnAdm(z)[j] * dH(x)[i]  + cnAdm(y)[j] * dH(y)[i]
                                                                                                                            - rnAdm(z)[j] * dH(z)[i];
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(z)[i], LIFT, fz);
            }
        }
    } else {
        // ---------- Penalized flux ------------------------------------------
        if (useResForUpw) {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = + nAdm(y)[j]  * dE(z)[i] - nAdm(z)[j]  * dE(y)[i]
                                                                           + (cnAdm(x)[j] * dresH(y)[i] + cnAdm(z)[j] * dresH(z)[i]
                                                                                                                                 - rnAdm(x)[j] * dresH(x)[i]) * upwinding;
                    fy[k] = + nAdm(z)[j]  * dE(x)[i] - nAdm(x)[j]  * dE(z)[i]
                                                                           + (cnAdm(y)[j] * dresH(z)[i] + cnAdm(x)[j] * dresH(x)[i]
                                                                                                                                 - rnAdm(y)[j] * dresH(y)[i]) * upwinding;
                    fz[k] = + nAdm(x)[j]  * dE(y)[i]  - nAdm(y)[j]  * dE(x)[i]
                                                                            + (cnAdm(z)[j] * dresH(x)[i]  + cnAdm(y)[j] * dresH(y)[i]
                                                                                                                                   - rnAdm(z)[j] * dresH(z)[i]) * upwinding;
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(z)[i], LIFT, fz);
            }
        } else {
#pragma omp parallel for private(i,j,k,e,f,fx,fy,fz)
            for (e = e1; e < e2; e++) {
                i = e * nfpfaces;
                for (k = 0; k < nfpfaces; k++) {
                    f = k / nfp;
                    j = e * faces + f;
                    fx[k] = + nAdm(y)[j]  * dE(z)[i] - nAdm(z)[j]  * dE(y)[i]
                                                                           + (cnAdm(x)[j] * dH(y)[i] + cnAdm(z)[j] * dH(z)[i]
                                                                                                                           - rnAdm(x)[j] * dH(x)[i]) * upwinding;
                    fy[k] = + nAdm(z)[j]  * dE(x)[i] - nAdm(x)[j]  * dE(z)[i]
                                                                           + (cnAdm(y)[j] * dH(z)[i] + cnAdm(x)[j] * dH(x)[i]
                                                                                                                           - rnAdm(y)[j] * dH(y)[i]) * upwinding;
                    fz[k] = + nAdm(x)[j]  * dE(y)[i]  - nAdm(y)[j]  * dE(x)[i]
                                                                            + (cnAdm(z)[j] * dH(x)[i]  + cnAdm(y)[j] * dH(y)[i]
                                                                                                                             - rnAdm(z)[j] * dH(z)[i]) * upwinding;
                    i++;
                }
                i = e * np;
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(x)[i], LIFT, fx);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(y)[i], LIFT, fy);
                add_m_v_prod<Real,np,nfpfaces>(&rhsH.set(z)[i], LIFT, fz);
            }
        }
    }
}
void DGExplicit::addCurvedFluxesToRHSElectric(
        const size_t e1,
        const size_t e2,
        const bool useResForUpw) {
    size_t c;
#pragma omp parallel for private(c)
    for (c = 0; c < nCurvedFaces; c++) {
        if (e1 <= curveFace[c].solverPosition
                && curveFace[c].solverPosition < e2) {
            curveFace[c].addFluxToRHSElectric(upwinding, useResForUpw);
        }
    }
}

void DGExplicit::addCurvedFluxesToRHSMagnetic(
        const size_t e1,
        const size_t e2,
        const bool useResForUpw) {
    size_t c;
#pragma omp parallel for private(c)
    for (c = 0; c < nCurvedFaces; c++) {
        if (e1 <= curveFace[c].solverPosition
                && curveFace[c].solverPosition < e2) {
            curveFace[c].addFluxToRHSMagnetic(upwinding, useResForUpw);
        }
    }
}

void DGExplicit::addRHSToFieldsElectric(
        const size_t e1,
        const size_t e2,
        const Real rkdt) {
    size_t init = getIndexOfElement(e1);
    size_t end = getIndexOfElement(e2);
    E.addProd_omp(init, end, getRHSElectric(), rkdt);
}

void DGExplicit::addRHSToFieldsMagnetic(
        const size_t e1,
        const size_t e2,
        const Real rkdt) {
    size_t init = getIndexOfElement(e1);
    size_t end = getIndexOfElement(e2);
    H.addProd_omp(init, end, getRHSMagnetic(), rkdt);
}

void DGExplicit::allocateRHSAndJumps() {
    rhsE.setSize(getFieldDOFs()/3);
    rhsH.setSize(getFieldDOFs()/3);
    dE.setSize(nK*nfp*4);
    dH.setSize(nK*nfp*4);
}

void DGExplicit::allocateMaps() {
    ExP = new Real**[nK];
    EyP = new Real**[nK];
    EzP = new Real**[nK];
    HxP = new Real**[nK];
    HyP = new Real**[nK];
    HzP = new Real**[nK];
    for (size_t e = 0; e < nK; e++) {
        ExP[e] = new Real*[faces];
        EyP[e] = new Real*[faces];
        EzP[e] = new Real*[faces];
        HxP[e] = new Real*[faces];
        HyP[e] = new Real*[faces];
        HzP[e] = new Real*[faces];
        for (size_t f = 0; f < faces; f++) {
            ExP[e][f] = NULL;
            EyP[e][f] = NULL;
            EzP[e][f] = NULL;
            HxP[e][f] = NULL;
            HyP[e][f] = NULL;
            HzP[e][f] = NULL;
        }
    }
    map_ = new Int**[nK];
    for (size_t e = 0; e < nK; e++) {
        map_[e] = new Int*[faces];
        for (size_t f = 0; f < faces; f++) {
            map_[e][f] = NULL;
        }
    }
}

void DGExplicit::assignMatrices(const CellGroup& cells) {
    // Assigns matrices.
    Cx = new const Real*[nK];
    Cy = new const Real*[nK];
    Cz = new const Real*[nK];
    size_t e;
#	ifdef SOLVER_DEDUPLICATE_OPERATORS
#	pragma omp parallel for private(e)
    for (e = 0; e < nK; e++) {
        array<StaMatrix<Real,np,np>,3> C;
        const ElemId id = cells.getIdOfRelPos(e);
        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
        C = cell->getCMatrices();
        for (size_t i = 0; i < 3; i++) {
            set<StaMatrix<Real,np,np>, lexCompareMat>::iterator it;
            it = CList.find(C[i]);
            if (it != CList.end()) {
                if (i == 0) {
                    Cx[e] = it->val();
                }
                if (i == 1) {
                    Cy[e] = it->val();
                }
                if (i == 2) {
                    Cz[e] = it->val();
                }
            } else {
                C[i].printInfo();
                throw Error("The following matrix is not stored.");
            }
        }
    }
    cout<< "Matrix deduplication: " << nK*3 << " --> " << CList.size() << endl;
#	else
    CList = new StaMatrix<Real,np,np>[3*nK];
#	pragma omp parallel for private(e)
    for (e = 0; e < nK; e++) {
        size_t id = cells.getIdOfRelPos(e);
        const CellTet<ORDER_N>* cell_ = cells.getPtrToCellWithId(id);
        for (size_t i = 0; i < 3; i++) {
            size_t j = 3 * e + i;
            if (i == 0) {
                Cx[e] = CList[j].val();
            }
            if (i == 1) {
                Cy[e] = CList[j].val();
            }
            if (i == 2) {
                Cz[e] = CList[j].val();
            }
        }
    }
#	endif
}

void DGExplicit::assignPointersToNeighbours(
        const CellGroup& cells,
        const Connectivities& map,
        const Mesh::Volume& mesh) {
    size_t nNeighs = 0;
    for (size_t k = 0; k < nK; k++) {
        const VolR* vol = cells(k)->getBase();
        for (size_t f = 0; f < faces; f++) {
            Face neigh = map.getNeighFace(Face(vol,f));
            ElemId id2;
            if (map.isDomainBoundary(neigh)) {
                id2 = vol->getId();
            } else {
                id2 = neigh.first->getId();
            }
            if (cells.isLocalId(id2)) {
                // Assigns ptrs to local cells and counts non local neigh.
                size_t e2 = cells.getRelPosOfId(id2);
                ExP[k][f] = &E.set(x)[e2 * np];
                EyP[k][f] = &E.set(y)[e2 * np];
                EzP[k][f] = &E.set(z)[e2 * np];
                HxP[k][f] = &H.set(x)[e2 * np];
                HyP[k][f] = &H.set(y)[e2 * np];
                HzP[k][f] = &H.set(z)[e2 * np];
            } else {
                nNeighs++;
            }
        }
    }
    // Stores neighbours information and allocates fields.
    vector<ElemId> neighId;
    neighId.reserve(nNeighs);
    nE.setSize(nNeighs*np);
    nH.setSize(nNeighs*np);
    nE.setAll((Real) 0.0);
    nH.setAll((Real) 0.0);
    size_t neigh = 0;
    for (size_t k = 0; k < nK; k++) {
        const VolR* vol = cells(k)->getBase();
        for (size_t f = 0; f < faces; f++) {
            Face neighF = map.getNeighFace(Face(vol,f));
            ElemId id2;
            if (map.isDomainBoundary(neighF)) {
                id2 = vol->getId();
            } else {
                id2 = neighF.first->getId();
            }
            if (!cells.isLocalId(id2)) {
                neighId.push_back(id2);
                ExP[k][f] = &nE.set(x)[neigh * np];
                EyP[k][f] = &nE.set(y)[neigh * np];
                EzP[k][f] = &nE.set(z)[neigh * np];
                HxP[k][f] = &nH.set(x)[neigh * np];
                HyP[k][f] = &nH.set(y)[neigh * np];
                HzP[k][f] = &nH.set(z)[neigh * np];
                neigh++;
            }
        }
    }
    //
    comm->initNeighbourFields(neighId);
    assert(neigh == nNeighs);
    assert(checkPtrsToNeigh());
}

vector<const BoundaryCondition*> DGExplicit::removeNonLocalBCs(
        const CellGroup* cells,
        const vector<const BoundaryCondition*>& bc) const {
    vector<const BoundaryCondition*> res;
    res.reserve(bc.size());
    for (uint i = 0; i < bc.size(); i++) {
        ElemId id = bc[i]->getCell()->getId();
        if (cells->isLocalId(id)) {
            res.push_back(bc[i]);
        }
    }
    return res;
}

void DGExplicit::BCToLocalArray(
        const BCGroup& bc,
        const CellGroup& cells,
        const Connectivities& map) {
    // ----------- SMA ------------------------------------------------
    // Counts SMAs and allocates, em boundaries are also considered.
    {
        vector<const BoundaryCondition*> smaPtr = bc.getPtrsToSMA();
        vector<const BoundaryCondition*> emPtr = bc.getPtrsToEMSourceBC();
        // Removes non local elements.
        smaPtr = removeNonLocalBCs(&cells, smaPtr);
        emPtr = removeNonLocalBCs(&cells, emPtr);
        // Stores em sources at boundaries.
        vector<const BoundaryCondition*> emAtDomainBound;
        for (size_t i = 0; i < emPtr.size(); i++) {
            if (map.isDomainBoundary(emPtr[i]->getCellFace())) {
                emAtDomainBound.push_back(emPtr[i]);
            }
        }
        nSMA = smaPtr.size() + emAtDomainBound.size();
        SMAe = new size_t[nSMA];
        SMAf = new size_t[nSMA];
        // Stores solver relative positions and faces.
        size_t j = 0;
        for (size_t i = 0; i < smaPtr.size(); i++) {
            SMAe[j] = cells.getRelPosOfId(smaPtr[i]->getCell()->getId());
            SMAf[j] = smaPtr[i]->getFace();
            j++;
        }
        for (size_t i = 0; i < emAtDomainBound.size(); i++) {
            SMAe[j] = cells.getRelPosOfId(emAtDomainBound[i]->getCell()->getId());
            SMAf[j] = emAtDomainBound[i]->getFace();
            j++;
        }
        assert(nSMA == j);
    }
    // ----------- PEC ------------------------------------------------
    // Counts and allocates.
    {
        vector<const BoundaryCondition*> pecPtr = bc.getPtrsToPEC();
        pecPtr = removeNonLocalBCs(&cells, pecPtr);
        nPEC = pecPtr.size();
        PECe = new size_t[nPEC];
        PECf = new size_t[nPEC];
        // Stores solver rel pos and faces.
        for (size_t i = 0; i < pecPtr.size(); i++) {
            PECe[i] = cells.getRelPosOfId(pecPtr[i]->getCell()->getId());
            PECf[i] = pecPtr[i]->getFace();
        }
    }
    // ----------- PMC ------------------------------------------------
    // Counts and allocates.
    {
        vector<const BoundaryCondition*> pmcPtr = bc.getPtrsToPMC();
        pmcPtr = removeNonLocalBCs(&cells, pmcPtr);
        nPMC = pmcPtr.size();
        PMCe = new size_t[nPMC];
        PMCf = new size_t[nPMC];
        // Stores solver rel pos and faces.
        for (size_t i = 0; i < pmcPtr.size(); i++) {
            PMCe[i] = cells.getRelPosOfId(pmcPtr[i]->getCell()->getId());
            PMCf[i] = pmcPtr[i]->getFace();
        }
    }
    //    {
    //        for (size_t i = 0; i < smb_->pMGroup->countSIBC(); i++) {
    //            const PMSurfaceSIBC* m =
    //                    dynamic_cast<const PMSurfaceSIBC*>(smb_->pMGroup->getPMSurface(i));
    //            if (m != NULL) {
    //                vector<const BoundaryCondition*> sibcPtr
    //                = bc.getMatId(m->getId());
    //                sibcPtr = removeNonLocalBCs(&cells, sibcPtr);
    //                dispersive.push_back(
    //                        new DGSIBC(*m, sibcPtr, cells, map, vmapM,
    //                                ExP, EyP, EzP, HxP, HyP, HzP));
    //            }
    //        }
    //    }
}

void DGExplicit::buildEMSources(
        const EMSourceGroup& em,
        const BCGroup& bc,
        const Connectivities& maps,
        const CellGroup& cells) {
    // Copies the sources structure into solver.
    for (size_t i = 0; i < em.getOf<PlaneWave>().size(); i++) {
        vector<const BoundaryCondition*> aux = bc.getPtrsToEMSourceBC();
        source.push_back(new DGPlaneWave(*(em(i)->castTo<PlaneWave>()), bc,
                maps, cells, comm, dE, dH, vmapM));
    }
    //    for (size_t i = 0; i < em.countDipoles(); i++) {
    //        vector<const BoundaryCondition*> aux = bc.get(Condition::emSource);
    //        source.push_back(
    //                new DGDipole(*em.getDipole(i), aux, maps, cells, dE, dH, vmapM));
    //    }
    //    for (size_t i = 0; i < em.countWaveports(); i++) {
    //        vector<const BoundaryCondition*> aux =	 bc.get(Condition::emSource);
    //        Waveport::Shape shape = em.getWaveport(i)->getShape();
    //        if (shape == Waveport::rectangular) {
    //            source.push_back(new DGWaveportRectangular(
    //                    *em.getWaveport(i), aux, maps, cells, dE, dH, vmapM));
    //        } else {
    //           throw Error("Unreckognized waveport shape.");
    //        }
    //    }
}



void DGExplicit::buildScalingFactors(
        const CellGroup& cells,
        const Connectivities& map) {
    buildFieldScalingFactors(cells);
    buildFluxScalingFactors(cells, map);
    buildCurvedFluxScalingFactors(cells, map);
}

void DGExplicit::buildCurvedFluxScalingFactors(
        const CellGroup& cells,
        const Connectivities& map) {
    // Counts curved faces.
    nCurvedFaces = 0;
    for (size_t e = 0; e < nK; e++) {
        ElemId id = cells.getIdOfRelPos(e);
        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
        for (size_t f = 0; f < cell->getFaces(); f++) {
            if (cell->isCurvedFace(f)) {
                nCurvedFaces++;
            }
        }
    }
    curveFace = new DGCurvedFace[nCurvedFaces];
    // Supress linear fluxes operators. Computes curved operators.
    size_t face = 0;
    for (size_t e = 0; e < nK; e++) {
        const VolR* vol = cells(e)->getBase();
        const CellTet<ORDER_N>* cell = cells.getPtrToCell(vol);
        Real impM, admM, impP, admP, impAv, admAv;
        for (size_t f = 0; f < faces; f++) {
            if (cell->isCurvedFace(f)) {
                // Builds CurvedFace information
                impM = cell->material->getImpedance();
                admM = cell->material->getAdmitance();
                ElemId nId = map.getNeighFace(Face(vol, f)).first->getId();
                const CellTet<ORDER_N>* neigh = cells.getPtrToCellWithId(nId);
                impP = neigh->material->getImpedance();
                admP = neigh->material->getAdmitance();
                impAv = (impM + impP) * 0.5;
                admAv = (admM + admP) * 0.5;
                curveFace[face++] = DGCurvedFace(
                        cell, f, e, rhsE, rhsH, dE, dH, dresE, dresH,
                        impP, admP, impAv, admAv);
                size_t i = e * faces + f;
                // Sets flux scaling factors to zero.
                CVecR3 zero(0.0, 0.0, 0.0);
                nAdm.set(i, zero);
                nImp.set(i, zero);
                cnAdm.set(i, zero);
                cnImp.set(i, zero);
                rnAdm.set(i, zero);
                rnImp.set(i, zero);
            }
        }
    }
}

void DGExplicit::buildMaterials(
        const CellGroup& cells,
        const OptionsSolverDGTD& arg) {
    //    // Creates Dispersive materials vars parameters and stores ptrs.
    //    const GroupPhysicalModels<PMVolumeDispersive> dispersives =
    //            smb_->pMGroup->getOf<PMVolumeDispersive>();
    //    for (size_t i = 0; i < dispersives.size(); i++) {
    //        dispersive.push_back(new DGDispersiveVolumic(*dispersives(i), cells));
    //    }
    //    // Creates PML materials variables parameters and stores pointers.
    //    const GroupPhysicalModels<PMVolumePML> pmls =
    //            smb_->pMGroup->getOf<PMVolumePML>();
    //    for (size_t i = 0; i < pmls.size(); i++) {
    //        const bool isConstCond = arg->isPMLConstantConductivityProfile();
    //        const Real cond = arg->getPMLConductivity();
    //        switch (pmls(i)->getOrientation()) {
    //        case PMVolumePML::Orientation::PMLx:
    //            dispersive.push_back(new DGPMLx(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        case PMVolumePML::Orientation::PMLy:
    //            dispersive.push_back(new DGPMLy(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        case PMVolumePML::Orientation::PMLz:
    //            dispersive.push_back(new DGPMLz(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        case PMVolumePML::Orientation::PMLxy:
    //            dispersive.push_back(new DGPMLxy(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        case PMVolumePML::Orientation::PMLyz:
    //            dispersive.push_back(new DGPMLyz(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        case PMVolumePML::Orientation::PMLzx:
    //            dispersive.push_back(new DGPMLzx(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        case PMVolumePML::Orientation::PMLxyz:
    //            dispersive.push_back(new DGPMLxyz(*pmls(i),cells,isConstCond,cond));
    //            break;
    //        default:
    //            break;
    //        }
    //    }
}

bool DGExplicit::checkPtrsToNeigh() const {
    bool res = true;
    for (size_t e = 0; e < nK; e++) {
        for (size_t f = 0; f < faces; f++) {
            bool problem = false;
            problem |= (ExP[e][f] == NULL);
            problem |= (EyP[e][f] == NULL);
            problem |= (EzP[e][f] == NULL);
            problem |= (HxP[e][f] == NULL);
            problem |= (HyP[e][f] == NULL);
            problem |= (HzP[e][f] == NULL);
            if (problem) {
                res = false;
                cout<< "ERROR@checkPtrsToNeigh" << endl;
                cout<< "Ptr to neigh " << e << "face " << f <<
                        " not assigned in solver init." << endl;
            }
        }
    }
    return res;
}

void DGExplicit::deduplicateVMaps(const CellGroup& cells) {
    // --- Copies vmapM -----------------------------------------------
    SimplexTet<ORDER_N> tet;
    for (size_t f = 0; f < faces; f++) {
        for (size_t i = 0; i < nfp; i++) {
            vmapM[f][i] = tet.sideNode(f,i);
        }
    }
    // --- Copies vmapP -----------------------------------------------
    // Allocates space for the deduplicated vmapP and inits all to -1.
    for (Int i = 0; i < 16; i++) {
        for (Int j = 0; j < Int(nfp); j++) {
            vmapP[i][j] = -1;
        }
    }
    // deduplicates vmapP in a list and points cell->vmap to them.
    for (size_t e = 0; e < nK; e++) {
        ElemId id = cells.getIdOfRelPos(e);
        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
        for (size_t f = 0; f < faces; f++) {
            // Checks if the vmapP[f] vector is in the list.
            bool found;
            for (Int i = 0; i < 16; i++) {
                found = true;
                for (size_t j = 0; j < nfp; j++)
                    if (vmapP[i][j] != Int(cell->vmapP[f][j])) {
                        found = false;
                        break;
                    }
                // Make the elements to point to their corresponding vmapP.
                if (found) {
                    map_[e][f] = vmapP[i];
                    break;
                }
            }
            // If not in the list, is added to the first empty row.
            if (!found) {
                for (Int i = 0; i < 16; i++) {
                    if (vmapP[i][0] == -1) {
                        for (size_t j = 0; j < nfp; j++) {
                            vmapP[i][j] = cell->vmapP[f][j];
                        }
                        // Points cell->vmapP to its corresponding one.
                        map_[e][f] = vmapP[i];
                        break;
                    }
                }
            }
        }
    }
    // Checks that all the elem vmapP pointers have been correctly addressed.
    bool problem = false;
    for (size_t e = 0; e < nK; e++) {
        for (size_t f = 0; f < faces; f++) {
            if (map_[e][f] == NULL) {
                problem = true;
                if (problem) {
                    cerr << endl << "ERROR: Solver::deduplicateVMaps()"
                            << "An vmapP has not been set. Elem "
                            << e << ", face " << f << endl;
                }
            }
        }
    }
    // Checks that vmapP has the correct values.
    for (size_t e = 0; e < nK; e++) {
        ElemId id = cells.getIdOfRelPos(e);
        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
        for (size_t f = 0; f < faces; f++) {
            for (size_t i = 0; i < nfp; i++) {
                if (Int(cell->vmapP[f][i]) != map_[e][f][i]) {
                    if (!problem) {
                        cerr << endl << "ERROR: Solver::deduplicateVMaps" << endl;
                        cerr << endl << "Check vmapsP" << endl;
                    }
                    cerr << endl << "Elem " << e << ", face " << f << endl;
                    problem = true;
                }
            }
        }
    }
}

//vector<const BoundaryCondition*>
//DGExplicit::removeNonLocalBCs(
//        const CellGroup* cells,
//        const vector<const BoundaryCondition*>& bc) const {
//    vector<const BoundaryCondition*> res;
//    res.reserve(bc.size());
//    for (size_t i = 0; i < bc.size(); i++) {
//        size_t id = bc[i]->getCell()->getId();
//        if (cells->isLocalId(id)) {
//            res.push_back(bc[i]);
//        }
//    }
//    return res;
//}

void DGExplicit::LTSSaveFieldsAndResidues(
        const size_t fKSave,
        const size_t lKSave) {
    const size_t init = fKSave * np;
    const size_t end = lKSave * np;
    savedE.copy(init, end, E);
    savedH.copy(init, end, H);
    savedResE.copy(init, end, resE);
    savedResH.copy(init, end, resH);
}

void DGExplicit::LTSLoadFieldsAndResidues(
        const size_t fKLoad,
        const size_t lKLoad) {
    const size_t init = fKLoad * np;
    const size_t end = lKLoad * np;
    E.copy(init, end, savedE);
    H.copy(init, end, savedH);
    resE.copy(init, end, savedResE);
    resH.copy(init, end, savedResH);
}


void DGExplicit::allocateFieldsForLTS() {
    dresE.setSize(nK*nfp*4);
    dresH.setSize(nK*nfp*4);
    savedE.setSize(getFieldDOFs()/3);
    savedH.setSize(getFieldDOFs()/3);
    savedResE.setSize(getFieldDOFs()/3);
    savedResH.setSize(getFieldDOFs()/3);
}

void DGExplicit::copyJumpsToResidueJumps(
        const size_t e1,
        const size_t e2) {
    size_t i, j, f, e;
#pragma omp parallel for private(e,i,f,j)
    for (e = e1; e < e2; e++) {
        i = e * nfp * faces;
        for (f = 0; f < faces; f++) {
            for (j = 0; j < nfp; j++) {
                dresE.set(x)[i] = dE(x)[i];
                dresE.set(y)[i] = dE(y)[i];
                dresE.set(z)[i] = dE(z)[i];
                dresH.set(x)[i] = dH(x)[i];
                dresH.set(y)[i] = dH(y)[i];
                dresH.set(z)[i] = dH(z)[i];
                i++;
            }
        }
    }
}


void DGExplicit::addRHSToResidueElectric(
        const size_t e1, const size_t e2,	const Real rkdt) {
    size_t i, j, e;
#pragma omp parallel for private(i,j,e)
    for (e = e1; e < e2; e++) {
        i = e * np;
        for (j = 0; j < np; j++) {
            resE.set(x)[i] += rhsE(x)[i] * rkdt;
            resE.set(y)[i] += rhsE(y)[i] * rkdt;
            resE.set(z)[i] += rhsE(z)[i] * rkdt;
            i++;
        }
    }
}

void DGExplicit::addRHSToResidueMagnetic(
        const size_t e1, const size_t e2, const Real rkdt) {
    size_t i, j, e;
#pragma omp parallel for private(i,j,e)
    for (e = e1; e < e2; e++) {
        i = e * np;
        for (j = 0; j < np; j++) {
            resH.set(x)[i] += rhsH(x)[i] * rkdt;
            resH.set(y)[i] += rhsH(y)[i] * rkdt;
            resH.set(z)[i] += rhsH(z)[i] * rkdt ;
            i++;
        }
    }
}

void DGExplicit::updateFieldsWithRes(
        const size_t e1,
        const size_t e2,
        const Real rkb) {
    updateFieldsWithResBase(e1,e2,rkb);
    for (size_t d = 0; d < dispersive.size(); d++) {
        dispersive[d]->updateWithRes(e1,e2,rkb);
    }
}

void DGExplicit::addRHSToRes(
        const size_t e1,
        const size_t e2,
        const Real rka,
        const Real dt) {
    size_t i = getIndexOfElement(e1);
    size_t j = getIndexOfElement(e2);
    resE.prod_omp(i,j, rka);
    resE.addProd_omp(i,j, rhsE, dt);
    resH.prod_omp(i,j, rka);
    resH.addProd_omp(i,j, rhsH, dt);
    // Polarization currents in dispersive materials.
    for (size_t d = 0; d < dispersive.size(); d++) {
        dispersive[d]->addRHSToRes(e1,e2,rka,dt);
    }
}
