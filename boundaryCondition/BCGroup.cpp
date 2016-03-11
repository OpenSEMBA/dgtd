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
 * BCGroup.cpp
 *
 *  Created on: Jul 8, 2013
 *      Author: luis
 */

#include "../../dgtd/core/BCGroup.h"

BCGroup::GroupBoundaryConditions(
        const Mesh::Volume& mesh,
        const EMSourceGroup& em,
        const PMGroup& pm,
        const CellGroup& cells,
        const Connectivities& map) {
    buildEMSourceBC(mesh, em, cells);
    buildPhysicalModelBC(mesh, pm, cells, map);
    removeOverlapped();
    check();
}

BCGroup& BCGroup::operator=(const BCGroup &rhs) {
    if (this == &rhs) {
        return *this;
    }
    embc = rhs.embc;
    pmbc = rhs.pmbc;
    sibc = rhs.sibc;
    return *this;
}

void BCGroup::buildEMSourceBC(
        const Mesh::Volume& mesh,
        const EMSourceGroup& em,
        const CellGroup& cells) {
    for (UInt i = 0; i < em.size(); i++) {
        // If emSource it has been defined using a layer it is transferred to
        // the border of elems within.
        vector<Face> border;
        if (em(i)->elems().size() == 1) {
            const ElemR* elem = em(i)->elems()(0)->castTo<ElemR>();
            if (!elem->is<HexR8>() || elem->getMatId() != MatId(0)) {
                throw Error("Invalid definition of Planewave.");
            }
            BoxR3 box = elem->getBound();
            GroupElements<const ElemR> elems = mesh.elems().getInsideBound(box);
            elems.removeMatId(MatId(0));
            border = mesh.getInternalBorder(elems.getOf<VolR>());
        } else {
            border = mesh.getInternalBorder(em(i)->elems());
        }

        // Builds boundary conditions.
        for (UInt j = 0; j < border.size(); j++) {
            const CellTet<ORDER_N>* auxCell = cells.getPtrToCell(border[j].first);
            UInt face = border[j].second;
            EMSourceBC auxBC(auxCell, face, em(i));
            embc.push_back(auxBC);
        }
    }
}

void BCGroup::buildPhysicalModelBC(
        const Mesh::Volume& mesh,
        const PMGroup& pm,
        const CellGroup& cells,
        const Connectivities& map) {
    Group<const SurfR> surf = mesh.elems().getOf<SurfR>();
    for (UInt i = 0; i < surf.size(); i++) {
        if (surf(i)->getMatId() !=  MatId(0)) {
            const PhysicalModel* mat = pm.getId(surf(i)->getMatId());
            Face tFace = map.getInnerFace(surf(i));
            if (tFace.first == NULL) {
                tFace = map.getOuterFace(surf(i));
            }
            if (tFace.first == NULL) {
                surf(i)->printInfo();
                throw Error("Surface with mat defined is floating.");
            }
            const CellTet<ORDER_N>* cell = cells.getPtrToCell(tFace.first);
            UInt face = tFace.second;
            if (!mat->is<PMSurfaceSIBC>()) {
                const PMPredefined* pred = mat->castTo<PMPredefined>();
                pmbc.push_back(PhysicalModelBC(cell, face, pred));
                if (!map.isDomainBoundary(tFace)) {
                    tFace = map.getOuterFace(surf(i));
                    cell = cells.getPtrToCell(tFace.first);
                    face = tFace.second;
                    pmbc.push_back(PhysicalModelBC(cell, face, pred));
                }
            } else {
                const PMSurfaceSIBC* matSibc = mat->castTo<PMSurfaceSIBC>();
                Face neigh = map.getNeighFace(tFace);
                const CellTet<ORDER_N>* nCell = cells.getPtrToCell(neigh.first);
                const UInt nFace = neigh.second;
                if (cell->isLocalSide(face, surf(i))) {
                    sibc.push_back(
                            SurfaceImpedanceBC(cell,face,nCell,nFace,matSibc));
                } else {
                    sibc.push_back(
                            SurfaceImpedanceBC(nCell,nFace,cell,face,matSibc));
                }
            }
        }
    }
}

void BCGroup::removeOverlapped() {
    // Hierarchy in boundary conditions.
    // PEC, PMC, SMA, SIBC > EM Source
    // Builds separated lists.
    vector<BoundaryCondition*> pec, pmc, sma, sibc, em;
    for (UInt i = 0; i < pmbc.size(); i++) {
        if(pmbc[i].getCondition()->is<PMSMA>()) {
            sma.push_back(&pmbc[i]);
        } else if (pmbc[i].getCondition()->is<PMPEC>()) {
            pec.push_back(&pmbc[i]);
        } else if (pmbc[i].getCondition()->is<PMPMC>()) {
            pmc.push_back(&pmbc[i]);
        } else if (pmbc[i].getCondition()->is<PMSurfaceSIBC>()) {
            sibc.push_back(&pmbc[i]);
        }
    }
    em.reserve(embc.size());
    for (UInt i = 0; i < embc.size(); i++) {
        em.push_back(&embc[i]);
    }
    // Removes bc overlapping PEC. After this, PEC is cleaned.
    pmc = removeCommons(pmc, pec);
    sma = removeCommons(sma, pec);
    em = removeCommons(em, pec);
    // Cleans PMC
    sma = removeCommons(sma, pmc);
    em = removeCommons(em, pmc);
    // Cleans SMA
    em = removeCommons(em, sma);
    // Cleans SIBC
    em = removeCommons(em, sibc);
    // Rebuilds lists.
    vector<PhysicalModelBC> auxPMBC;
    auxPMBC.reserve(pec.size() + pmc.size() + sma.size() + sibc.size());
    for (UInt i = 0; i < pec.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(pec[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    for (UInt i = 0; i < pmc.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(pmc[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    for (UInt i = 0; i < sma.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(sma[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    for (UInt i = 0; i < sibc.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(sibc[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    pmbc = auxPMBC;
    vector<EMSourceBC> auxEMBC;
    auxEMBC.reserve(em.size());
    for (UInt i = 0; i < em.size(); i++) {
        const EMSourceBC* aux = dynamic_cast<EMSourceBC*>(em[i]);
        auxEMBC.push_back(EMSourceBC(*aux));
    }
    embc = auxEMBC;
}

vector<BoundaryCondition*> BCGroup::removeCommons(
        const vector<BoundaryCondition*>& low,
        const vector<BoundaryCondition*>& high) const {
    vector<BoundaryCondition*> res;
    res.reserve(low.size());
    for (UInt i = 0; i < low.size(); i++) {
        bool isPresentInHigh = false;
        for (UInt j = 0; j < high.size(); j++) {
            if (low[i]->hasSameBoundary(*high[j])) {
                isPresentInHigh = true;
                break;
            }
        }
        if (!isPresentInHigh) {
            res.push_back(low[i]);
        }
    }
    return res;
}

vector<const BoundaryCondition*> BCGroup::getPtrsToPEC() const {
    vector<const BoundaryCondition*> res;
    for (UInt i = 0; i < pmbc.size(); i++) {
        if (pmbc[i].getCondition()->is<PMPEC>()) {
            const BoundaryCondition* ptr = &pmbc[i];
            res.push_back(ptr);
        }
    }
    return res;
}

vector<const BoundaryCondition*> BCGroup::getPtrsToPMC() const {
    vector<const BoundaryCondition*> res;
    for (UInt i = 0; i < pmbc.size(); i++) {
        if (pmbc[i].getCondition()->is<PMPMC>()) {
            const BoundaryCondition* ptr = &pmbc[i];
            res.push_back(ptr);
        }
    }
    return res;
}

vector<const BoundaryCondition*> BCGroup::getPtrsToSMA() const {
    vector<const BoundaryCondition*> res;
    for (UInt i = 0; i < pmbc.size(); i++) {
        if (pmbc[i].getCondition()->is<PMSMA>()) {
            const BoundaryCondition* ptr = &pmbc[i];
            res.push_back(ptr);
        }
    }
    return res;
}

vector<const BoundaryCondition*> BCGroup::getPtrsToSIBC() const {
    vector<const BoundaryCondition*> res;
    for (UInt i = 0; i < sibc.size(); i++) {
        const BoundaryCondition* ptr = &sibc[i];
        res.push_back(ptr);
    }
    return res;
}

vector<const BoundaryCondition*> BCGroup::getPtrsToEMSourceBC() const {
    vector<const BoundaryCondition*> res;
    for (UInt i = 0; i < embc.size(); i++) {
        const BoundaryCondition* ptr = &embc[i];
        res.push_back(ptr);
    }
    return res;
}

vector<const BoundaryCondition*> BCGroup::getPtrsToBC(const EMSourceBase* pw) const {
    vector<const BoundaryCondition*> res;
    res.reserve(embc.size());
    for (UInt i = 0; i < embc.size(); i++) {
        if (pw == embc[i].getCondition()) {
            res.push_back(&embc[i]);
        }
    }
    return res;
}

void BCGroup::printInfo() const {
    cout << "--- BCGroup info ---" << endl;
    cout << "EM Source BCs: " << embc.size() << endl;
    cout << "Physical Model BCs: " << pmbc.size() << endl;
    cout << "Surface Impedance BCs: " << sibc.size() << endl;
}

void BCGroup::check() const {
    checkEMSourcesAreSetInVacuum();
    //assert(checkOverlapping());
}

vector<const BoundaryCondition*> BCGroup::getPtrsToBCWithMatId(
        const MatId id) const {
    vector<const BoundaryCondition*> res;
    res.reserve(pmbc.size());
    for (UInt i = 0; i < pmbc.size(); i++) {
        if (pmbc[i].getCondition()->getId() == id) {
            res.push_back(&pmbc[i]);
        }
    }
    res.reserve(sibc.size());
    for (UInt i = 0; i < sibc.size(); i++) {
        if (sibc[i].getCondition()->getId() == id) {
            res.push_back(&sibc[i]);
        }
    }
    return res;
}

void BCGroup::checkEMSourcesAreSetInVacuum() const {
    UInt nBC = embc.size();
    for (UInt i = 0; i < nBC; i++) {
        if (!embc[i].getCell()->material->isVacuum()) {
            throw Error("ElectromagneticSource BC has been defined over a not vacuum cell.");
        }
    }
}

bool BCGroup::checkOverlapping() const {
    bool repeated = false;
    // Check repeated elements in embc and pmbc
    for (UInt i = 0; i < embc.size(); i++) {
        for (UInt j = 0; j < pmbc.size(); j++) {
            if (embc[i].getCell() == pmbc[j].getCell()
                    && embc[i].getFace()  == pmbc[j].getFace()) {
                repeated = true;
                break;
            }
        }
    }
    if (repeated) {
        throw Error("Overlapping boundary conditions detected.");
    }
    return repeated;
}
