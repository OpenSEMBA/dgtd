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
 * DGSource.cpp
 *
 *  Created on: Sep 2, 2013
 *      Author: luis
 */

#include "../../dg/sources/DGSource.h"

DGSource::DGSource() {
    dExT = NULL;
    dEyT = NULL;
    dEzT = NULL;
    dHxT = NULL;
    dHyT = NULL;
    dHzT = NULL;
    dExS = NULL;
    dEyS = NULL;
    dEzS = NULL;
    dHxS = NULL;
    dHyS = NULL;
    dHzS = NULL;
    dExTNB = NULL;
    dEyTNB = NULL;
    dEzTNB = NULL;
    dHxTNB = NULL;
    dHyTNB = NULL;
    dHzTNB = NULL;
}

DGSource::~DGSource() {

}

void DGSource::initSource(
        const BCGroup& bc,
        const Connectivities& map,
        const CellGroup& cells,
        FieldR3& dE,
        FieldR3& dH,
        const Int vmapM[faces][nfp]) {
    vector<pair<size_t, size_t> > total, scatt, totalNotBacked;
    total = getTotalFieldElemFaces(bc, map, cells);
    scatt = getScattFieldElemFaces(bc, map, cells);
    totalNotBacked = getTotalNotBackedFieldElemFaces(bc, map, cells);
    // Set fields to zero.
    ETInc.setAll((Real) 0.0);
    HTInc.setAll((Real) 0.0);
    ESInc.setAll((Real) 0.0);
    HSInc.setAll((Real) 0.0);
    EIncNB.setAll((Real) 0.0);
    HIncNB.setAll((Real) 0.0);
    const size_t nETF = total.size();
    const size_t nESF = scatt.size();
    const size_t nETFNB = totalNotBacked.size();
    // Allocates and sets jumps pointers.
    // The pointers point to the beginning of the face that they have to
    // update on each iteration.
    dExT = new Real*[nETF];
    dEyT = new Real*[nETF];
    dEzT = new Real*[nETF];
    dHxT = new Real*[nETF];
    dHyT = new Real*[nETF];
    dHzT = new Real*[nETF];
    for (size_t j = 0; j < nETF; j++) {
        size_t e = total[j].first;
        size_t f = total[j].second;
        size_t pos = e * nfp * faces + f * nfp;
        dExT[j] = &dE.set(x)[pos];
        dEyT[j] = &dE.set(y)[pos];
        dEzT[j] = &dE.set(z)[pos];
        dHxT[j] = &dH.set(x)[pos];
        dHyT[j] = &dH.set(y)[pos];
        dHzT[j] = &dH.set(z)[pos];
    }
    dExS = new Real*[nESF];
    dEyS = new Real*[nESF];
    dEzS = new Real*[nESF];
    dHxS = new Real*[nESF];
    dHyS = new Real*[nESF];
    dHzS = new Real*[nESF];
    for (size_t j = 0; j < nESF; j++) {
        size_t e = scatt[j].first;
        size_t f = scatt[j].second;
        size_t pos = e * nfp * faces + f * nfp;
        dExS[j] = &dE.set(x)[pos];
        dEyS[j] = &dE.set(y)[pos];
        dEzS[j] = &dE.set(z)[pos];
        dHxS[j] = &dH.set(x)[pos];
        dHyS[j] = &dH.set(y)[pos];
        dHzS[j] = &dH.set(z)[pos];
    }
    dExTNB = new Real*[nETFNB];
    dEyTNB = new Real*[nETFNB];
    dEzTNB = new Real*[nETFNB];
    dHxTNB = new Real*[nETFNB];
    dHyTNB = new Real*[nETFNB];
    dHzTNB = new Real*[nETFNB];
    for (size_t j = 0; j < nETFNB; j++) {
        size_t e = totalNotBacked[j].first;
        size_t f = totalNotBacked[j].second;
        size_t pos = e * (nfp*faces) + f * nfp;
        dExTNB[j] = &dE.set(x)[pos];
        dEyTNB[j] = &dE.set(y)[pos];
        dEzTNB[j] = &dE.set(z)[pos];
        dHxTNB[j] = &dH.set(x)[pos];
        dHyTNB[j] = &dH.set(y)[pos];
        dHzTNB[j] = &dH.set(z)[pos];
    }
    // List of elements.
    ETFe.resize(nETF);
    for (size_t i = 0; i < nETF; i++) {
        ETFe[i] = total[i].first;
    }
    ESFe.resize(nESF);
    for (size_t i = 0; i < nESF; i++) {
        ESFe[i] = scatt[i].first;
    }
    ETFNBe.resize(nETFNB);
    for (size_t i = 0; i < nETFNB; i++) {
        ETFNBe[i] = totalNotBacked[i].first;
    }
}

void DGSource::addJumps(
        const size_t e1,
        const size_t e2) {
    size_t j, k, pos;
    // Total field jumps.
    for (j = 0; j < ETInc.size(); j++) {
        if (e1 <= ETFe[j] && ETFe[j] < e2) {
            for (k = 0; k < nfp; k++) {
                pos = j * nfp + k;
                dExT[j][k] -= ETInc(x)[pos];
                dEyT[j][k] -= ETInc(y)[pos];
                dEzT[j][k] -= ETInc(z)[pos];
                dHxT[j][k] -= HTInc(x)[pos];
                dHyT[j][k] -= HTInc(y)[pos];
                dHzT[j][k] -= HTInc(z)[pos];
            }
        }
    }
    // Scatt field jumps.
    for (j = 0; j < ESInc.size(); j++) {
        if (e1 <= ESFe[j] && ESFe[j] < e2) {
            for (k = 0; k < nfp; k++) {
                pos = j * nfp + k;
                dExS[j][k] += ESInc(x)[pos];
                dEyS[j][k] += ESInc(y)[pos];
                dEzS[j][k] += ESInc(z)[pos];
                dHxS[j][k] += HSInc(x)[pos];
                dHyS[j][k] += HSInc(y)[pos];
                dHzS[j][k] += HSInc(z)[pos];
            }
        }
    }
    // Computes TFNB excitation jumps.
    for (j = 0; j < EIncNB.size(); j++) {
        if (e1 <= ETFNBe[j] && ETFNBe[j] < e2) {
            for (k = 0; k < nfp; k++) {
                pos = j * nfp + k;
                // Total field jumps of not backed elements.
                // The inc. field is substracted to computed by the SMA
                dExTNB[j][k] -= EIncNB(x)[pos];
                dEyTNB[j][k] -= EIncNB(y)[pos];
                dEzTNB[j][k] -= EIncNB(z)[pos];
                dHxTNB[j][k] -= HIncNB(x)[pos];
                dHyTNB[j][k] -= HIncNB(y)[pos];
                dHzTNB[j][k] -= HIncNB(z)[pos];
            }
        }
    }
}

CVecR3* DGSource::initPositions(
        const vector<pair<size_t, size_t> >& elemFace,
        const CellGroup& cells) const {
    const size_t nE = elemFace.size();
    CVecR3 *pos;
    pos = new CVecR3 [nE * nfp];
    for (size_t i = 0; i < nE; i++) {
        ElemId id = cells.getIdOfRelPos(elemFace[i].first);
        size_t f = elemFace[i].second;
        for (size_t j = 0; j < nfp; j++) {
            pos[i*nfp+j] =
                    cells.getPtrToCellWithId(id)->getSideNodePos(f,j);
        }
    }
    return pos;
}

vector<pair<size_t, size_t>> DGSource::getTotalFieldElemFaces(
        const BCGroup& bc,
        const Connectivities& map,
        const CellGroup& cells) const {
    vector<pair<size_t, size_t> > res;
    for (size_t i = 0; i < bc.embc.size(); i++) {
        if (!map.isDomainBoundary(bc.embc[i].getCellFace())) {
            const ElemId id1 = bc.embc[i].getCell()->getId();
            const size_t f1 = bc.embc[i].getFace();
            if (cells.isLocalId(id1)) {
                pair<size_t,size_t> aux1(cells.getRelPosOfId(id1), f1);
                res.push_back(aux1);
            }
        }
    }
    return res;
}

vector<pair<size_t, size_t>> DGSource::getScattFieldElemFaces(
        const BCGroup& bc,
        const Connectivities& map,
        const CellGroup& cells) const {
    vector<pair<size_t,size_t> > res;
    for (size_t i = 0; i <bc.embc.size(); i++) {
        Face bcFace = bc.embc[i].getCellFace();
        if (!map.isDomainBoundary(bcFace)) {
            Face outer = map.getNeighFace(bcFace);
            ElemId id2 = outer.first->getId();
            size_t f2 = outer.second;
            if (cells.isLocalId(id2)) {
                size_t e2 = cells.getRelPosOfId(id2);
                pair<size_t,size_t> aux2(e2, f2);
                res.push_back(aux2);
            }
        }
    }
    return res;
}

vector<pair<size_t, size_t>> DGSource::getTotalNotBackedFieldElemFaces(

        const BCGroup& bc,
        const Connectivities& map,
        const CellGroup& cells) const {
    vector<pair<size_t,size_t> > res;
    for (size_t i = 0; i < bc.embc.size(); i++) {
        Face inner = bc.embc[i].getCellFace();
        if (map.isDomainBoundary(inner)) {
            const ElemId id1 = inner.first->getId();
            const size_t f1 = inner.second;
            if (cells.isLocalId(id1)) {
                pair<size_t,size_t> aux1(cells.getRelPosOfId(id1), f1);
                res.push_back(aux1);
            }
        }
    }
    return res;
}

