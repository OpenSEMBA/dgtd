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
 * DGPlaneWave.cpp
 *
 *  Created on: Sep 11, 2012
 *      Author: luis
 */
#include "../../dg/sources/DGPlaneWave.h"

DGPlaneWave::DGPlaneWave(
        const PlaneWave& pw,
        const BCGroup& bc,
        const Connectivities& map,
        const CellGroup& cells,
        const Comm* comm,
        FieldR3& dE, FieldR3& dH,
        const Int vmapM[faces][nfp]) :
        PlaneWave(pw) {
    	initSource(bc, map, cells, dE, dH, vmapM);
    	initWaveNumberPosition(bc, map, cells, comm, vmapM);
}


DGPlaneWave::~DGPlaneWave() {
}

void DGPlaneWave::printInfo() const {
    cout << " --- DGPlaneWave Info ---" << endl;
    cout << "# ETF: " << ETInc.size() << endl;
    cout << "# ESF: " << ESInc.size() << endl;
    cout << "# ETFNB: " << EIncNB.size() << endl;
}

void DGPlaneWave::computeExcitation(
        const Real intTime, const Real minDT) {
    computeExcitationField(ETInc, HTInc, kNPosTF, ETInc.size(), intTime);
    computeExcitationField(ESInc, HSInc, kNPosSF, ESInc.size(), intTime);
    computeExcitationField(EIncNB, HIncNB, kNPosTFNB, EIncNB.size(), intTime);
}

void DGPlaneWave::computeExcitationField(
        FieldR3& EInc,
        FieldR3& HInc,
        const Real* vPos,
        const UInt nE,
        const Real intTime) {
    // Computes the plane wve excitation corresponding to the face f.
    // The face contins nfp nodes and it is computed for time intTime.
    const UInt nFields = nfp*nE;
    for (UInt j = 0; j < nFields; j++) {
        Real delayedTime = intTime - vPos[j];
        if (delayedTime >= 0) {
            pair<CVecR3,CVecR3> EHInc = getElectromagneticField(delayedTime);
            EInc(x)[j] = EHInc.first(0);
            EInc(y)[j] = EHInc.first(1);
            EInc(z)[j] = EHInc.first(2);
            HInc(x)[j] = EHInc.second(0);
            HInc(y)[j] = EHInc.second(1);
            HInc(z)[j] = EHInc.second(2);
        } else {
            EInc(x)[j] = (Real) 0.0;
            EInc(y)[j] = (Real) 0.0;
            EInc(z)[j] = (Real) 0.0;
            HInc(x)[j] = (Real) 0.0;
            HInc(y)[j] = (Real) 0.0;
            HInc(z)[j] = (Real) 0.0;
        }
    }
}

void DGPlaneWave::initWaveNumberPosition(
        const BCGroup& bc,
        const Connectivities& map,
        const CellGroup& cells,
        const Comm* comm,
        const Int vmapM[faces][nfp]) {
    // Generates kNPosTSF for PW excitation.
    // Computes krmin position.
    Real krmin = 0.0;
    bool krminSet = false;
    CVecR3 nodePos;
    // Total field.
    vector<pair<UInt, UInt> > total;
    total = getTotalFieldElemFaces(bc, map, cells);
    kNPosTF = new Real[ETInc.size() * nfp];
    for (UInt j = 0; j < ETInc.size(); j++) {
        ElemId id = cells.getIdOfRelPos(total[j].first);
        UInt f = total[j].second;
        UInt pos = j * nfp;
        for (UInt k = 0; k < nfp; k++) {
            const CellTet<N>* cell = cells.getPtrToCellWithId(id);
            nodePos = cell->n[vmapM[f][k]];
            kNPosTF[pos + k] = getWaveDirection().dot(nodePos) / Constants::c0;
            // Stores minimum kNPos value.
            if (!krminSet) {
                krminSet = true;
                krmin = kNPosTF[pos + k];
            }
            else if (kNPosTF[pos + k] < krmin)
                krmin = kNPosTF[pos + k];
        }
    }
    // Scattered field.
    vector<pair<UInt, UInt> > scatt;
    scatt = getScattFieldElemFaces(bc, map, cells);
    kNPosSF = new Real[ESInc.size() * nfp];
    for (UInt j = 0; j < ESInc.size(); j++) {
        ElemId id = cells.getIdOfRelPos(scatt[j].first);
        UInt f = scatt[j].second;
        UInt pos = j * nfp;
        for (UInt k = 0; k < nfp; k++) {
            const CellTet<N>* cell = cells.getPtrToCellWithId(id);
            nodePos = cell->n[vmapM[f][k]];
            kNPosSF[pos + k] = getWaveDirection().dot(nodePos) / Constants::c0;
            // Stores minimum kNPos value.
            if (!krminSet) {
                krminSet = true;
                krmin = kNPosSF[pos + k];
            }
            else if (kNPosSF[pos + k] < krmin) {
                krmin = kNPosSF[pos + k];
            }
        }
    }
    // Total field not backed.
    vector<pair<UInt, UInt> > totalNotBacked;
    totalNotBacked = getTotalNotBackedFieldElemFaces(bc, map, cells);
    kNPosTFNB = new Real[EIncNB.size() * nfp];
    for (UInt j = 0; j < EIncNB.size(); j++) {
        ElemId id = cells.getIdOfRelPos(totalNotBacked[j].first);
        UInt f = totalNotBacked[j].second;
        UInt pos = j * nfp;
        for (UInt k = 0; k < nfp; k++) {
            const CellTet<N>* cell = cells.getPtrToCellWithId(id);
            nodePos = cell->n[vmapM[f][k]];
            kNPosTFNB[pos + k] = getWaveDirection().dot(nodePos) / Constants::c0;
            // Stores minimum kNPos value.
            if (!krminSet) {
                krminSet = true;
                krmin = kNPosTFNB[pos + k];
            }
            else if (kNPosTFNB[pos + k] < krmin) {
                krmin = kNPosTFNB[pos + k];
            }
        }
    }
    // Syncs minimum.
    krmin = comm->reduceToGlobalMinimum(krmin);
    // Adds krmin to kNPosTSF and kNPosTFNB.
    for (UInt j = 0; j < ETInc.size() * nfp; j++) {
        kNPosTF[j] -= krmin;
    }
    for (UInt j = 0; j < ESInc.size() * nfp; j++) {
        kNPosSF[j] -= krmin;
    }
    for (UInt j = 0; j < EIncNB.size() * nfp; j++) {
        kNPosTFNB[j] -= krmin;
    }
}
