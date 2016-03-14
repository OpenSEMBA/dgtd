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
 * Group.cpp
 *
 *  Created on: Jul 8, 2013
 *      Author: luis
 */

#include "Group.h"

namespace SEMBA {
namespace Cudg3d {
namespace BoundaryCondition {

Group::Group(
        const Mesh::Volume& mesh,
        const SourceGroup& em,
        const PMGroup& pm) {
    buildEMSourceBC_(mesh, em);
    buildPhysicalModelBC_(mesh, pm);
    removeOverlapped();
    check();
}

void Group::buildEMSourceBC_(
        const Mesh::Volume& mesh,
        const SourceGroup& em) {
    for (size_t i = 0; i < em.size(); i++) {
        // If emSource it has been defined using a layer it is transferred to
        // the border of elems within.
        vector<Geometry::Element::Face> border;
        if (em(i)->elems().size() == 1) {
            const Geometry::ElemR* elem =
                    em(i)->elems()(0)->castTo<Geometry::ElemR>();
            if (!elem->is<Geometry::HexR8>() || elem->getMatId() != MatId(0)) {
                throw logic_error("Invalid definition of Planewave.");
            }
            Geometry::BoxR3 box = elem->getBound();
            Geometry::ElemRGroup elems = mesh.elems().getInsideBound(box);
            elems.removeMatId(MatId(0));
            border = mesh.getInternalBorder(elems.getOf<Geometry::VolR>());
        } else {
            border = mesh.getInternalBorder(em(i)->elems());
        }

        // Builds boundary conditions.
        for (size_t j = 0; j < border.size(); j++) {
            Geometry::Element::Face localFace = border[j];
            Geometry::Element::Face neighFace(NULL, (size_t) 0);
            add(new EMSourceBC(localFace, neighFace, em(i)));
        }
    }
}

void Group::buildPhysicalModelBC_(
        const Mesh::Volume& mesh,
        const PMGroup& pm) {
    Geometry::SurfRGroup surf = mesh.elems().getOf<Geometry::SurfR>();
    for (size_t i = 0; i < surf.size(); i++) {
        if (surf(i)->getModel() != NULL) {
            Geometry::Element::Face localFace =
                    mesh.getConnectivities()->getInnerFace(surf(i));
            Geometry::Element::Face neighFace =
                    mesh.getConnectivities()->getOuterFace(surf(i));
            if (localFace.first == NULL) {
                localFace = neighFace;
                neighFace = Geometry::Element::Face(NULL, (size_t) 0);
            }
            if (localFace.first == NULL) {
                surf(i)->printInfo();
                throw logic_error("Surface with mat defined is floating.");
            }
            const PhysicalModel::PhysicalModel* mat = surf(i)->getModel());
            add(new PhysicalModelBC(localFace, neighFace, mat));
        }
    }
}

void Group::removeOverlapped() {
    // Hierarchy in boundary conditions: PEC, PMC, SMA, SIBC > EM Source
    vector<BoundaryCondition*> pec, pmc, sma, sibc, em;
    for (size_t i = 0; i < pmbc.size(); i++) {
        if(pmbc[i].getCondition()->is<PhysicalModel::Predefined::SMA>()) {
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
    for (size_t i = 0; i < embc.size(); i++) {
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
    for (size_t i = 0; i < pec.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(pec[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    for (size_t i = 0; i < pmc.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(pmc[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    for (size_t i = 0; i < sma.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(sma[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    for (size_t i = 0; i < sibc.size(); i++) {
        const PhysicalModelBC* aux = dynamic_cast<PhysicalModelBC*>(sibc[i]);
        auxPMBC.push_back(PhysicalModelBC(*aux));
    }
    pmbc = auxPMBC;
    vector<EMSourceBC> auxEMBC;
    auxEMBC.reserve(em.size());
    for (size_t i = 0; i < em.size(); i++) {
        const EMSourceBC* aux = dynamic_cast<EMSourceBC*>(em[i]);
        auxEMBC.push_back(EMSourceBC(*aux));
    }
    embc = auxEMBC;
}

vector<Base*> Group::removeCommons(
        const vector<Base*>& low,
        const vector<Base*>& high) const {
    vector<BoundaryCondition::Base*> res;
    res.reserve(low.size());
    for (size_t i = 0; i < low.size(); i++) {
        bool isPresentInHigh = false;
        for (size_t j = 0; j < high.size(); j++) {
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

void Group::printInfo() const {
    cout << "--- Group info ---" << endl;
    cout << "EM Source BCs: " << embc.size() << endl;
    cout << "Physical Model BCs: " << pmbc.size() << endl;
    cout << "Surface Impedance BCs: " << sibc.size() << endl;
}

}
}
}
