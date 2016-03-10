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
 * CellGroup.cpp
 *
 *  Created on: Aug 29, 2012
 *      Author: luis
 */
#include "../../dgtd/core/CellGroup.h"

CellGroup::CellGroup(const Mesh::Volume& mesh, const PMGroup& pMGroup) {
    Group<const Tet> tet = mesh.elems().getOf<Tet>();
    cell.resize(tet.size(), NULL);
    cellOffsetId = tet(0)->getId().toUInt();
    // Reserves space for cell vectors.
    vector<const Tet*> linear, quadratic;
    for (UInt k = 0; k < tet.size(); k++) {
        if (!tet(k)->isCurved()) {
            linear.push_back(tet(k));
        } else {
            quadratic.push_back(tet(k));
        }
    }
    linTet.resize(linear.size(), CellTet4<ORDER_N>());
    quadTet.resize(quadratic.size(), CellTet10<ORDER_N>());
    for (UInt k = 0; k < linear.size(); k++) {
        linTet[k] = CellTet4<ORDER_N>(linear[k], pMGroup);
        cell[linTet[k].getId().toUInt() - cellOffsetId] = &linTet[k];
    }
    for (UInt k = 0; k < quadratic.size(); k++) {
        quadTet[k] = CellTet10<ORDER_N>(quadratic[k], pMGroup);
        cell[quadTet[k].getId().toUInt() - cellOffsetId] = &quadTet[k];
    }

    Connectivities map(mesh.elems());
    buildNodalMaps(map);
    check(map);
}

CellGroup::~CellGroup() {
    // TODO Auto-generated destructor stub
}

const CellTet<ORDER_N>* CellGroup::operator()(const UInt i) const {
    return cell[i];
}

const CellTet<ORDER_N>* CellGroup::getPtrToCell(const VolR* elem) const {
    return getPtrToCellWithId(elem->getId());
}

const CellTet<ORDER_N>* CellGroup::getPtrToCellWithId(const ElemId& id) const {
    return cell[id.toUInt() - cellOffsetId];
}

void CellGroup::buildNodalMaps(const Connectivities& map) {
    for (UInt e = 0; e < cell.size(); e++) {
        for (UInt f = 0; f < cell[e]->getFaces(); f++) {
            Face local(cell[e]->getBase(),f);
            Face neigh = map.getNeighFace(local);
            if (map.isDomainBoundary(local)) {
                neigh = local;
            }
            const CellTet<ORDER_N>* c2 = getPtrToCell(neigh.first);
            for (UInt i = 0; i < cell[e]->getNfp(); i++) {
                const CVecR3 posM = cell[e]->getSideNodePos(f,i);
                for (UInt j = 0; j < c2->getNfp(); j++) {
                    const UInt f2 = neigh.second;
                    const CVecR3 posP = c2->getSideNodePos(f2,j);
                    if (posM == posP) {
                        cell[e]->vmapP[f][i] = cell[e]->getSideNode(f2,j);
                        break;
                    }
                }
            }
        }
    }
}

void CellGroup::check(const Connectivities& map) const {
    checkNodalMaps(map);
}

void CellGroup::checkNodalMaps(const Connectivities& map) const {
    // Checks for vmap.
    UInt nK = cell.size();
    for (UInt e = 0; e < nK; e++) {
        for (int f = 0; f < 4; f++) {
            Face local(cell[e]->getBase(),f);
            Face neigh = map.getNeighFace(local);
            if (map.isDomainBoundary(local)) {
                neigh = local;
            }
            const CellTet<ORDER_N>* c2 = getPtrToCell(neigh.first);
            for (UInt i = 0; i < cell[e]->getNfp(); i++) {
                int neighNode = cell[e]->vmapP[f][i];
                if (cell[e]->getSideNodePos(f, i) != c2->n[neighNode]) {
                    cerr << "Elem " << e << ", face " << f << endl;
                    throw Error("vmapP contains errors.");
                }
            }
        }
    }
}

