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

#include "../../dgtd/core/BoundaryCondition.h"

BoundaryCondition::BoundaryCondition() {
    cell_ = NULL;
    face_ = 0;
}


BoundaryCondition::BoundaryCondition(const CellTet<ORDER_N>* cell, UInt face) {
    cell_ = cell;
    face_ = face;
}

BoundaryCondition::~BoundaryCondition() {

}

bool BoundaryCondition::hasSameBoundary(
        const BoundaryCondition& other) const {
    return (cell_ == other.cell_ && face_ == other.face_);
}

BoundaryCondition& BoundaryCondition::operator =(const BoundaryCondition& rhs) {
    if (this == &rhs) {
        return *this;
    }
    cell_ = rhs.cell_;
    face_ = rhs.face_;
    return *this;
}

Face BoundaryCondition::getCellFace() const {
    return Face(cell_->getBase(), face_);
}

void BoundaryCondition::printInfo() const {
    cout << "--- BC info ---" << endl;
    cout << "Cell Id:" << cell_->getId() << " Face: " << face_ << endl;
}

EMSourceBC::EMSourceBC() {
    em_ = NULL;
}

EMSourceBC::~EMSourceBC() {

}

EMSourceBC::EMSourceBC(
        const CellTet<ORDER_N>* e,
        const UInt f,
        const EMSourceBase* bc) :
            BoundaryCondition(e, f) {
    em_ = bc;
}

EMSourceBC& EMSourceBC::operator=(const EMSourceBC& rhs) {
    if (this == &rhs) {
        return *this;
    }
    BoundaryCondition::operator=(rhs);
    em_ = rhs.em_;
    return *this;
}

PhysicalModelBC::PhysicalModelBC() {
    pm_ = NULL;
}

PhysicalModelBC::~PhysicalModelBC() {

}

PhysicalModelBC::PhysicalModelBC(
        const CellTet<ORDER_N>* cell__,
        UInt face__,
        const PMPredefined* bc) :
            BoundaryCondition(cell__, face__) {
    pm_ = bc;
}

PhysicalModelBC& PhysicalModelBC::operator=(const PhysicalModelBC& rhs) {
    if (this == &rhs) {
        return *this;
    }
    BoundaryCondition::operator =(rhs);
    pm_ = rhs.pm_;
    return *this;
}

SurfaceImpedanceBC::SurfaceImpedanceBC() {
    cellD_ = NULL;
    faceD_ = 0;
    sibc_ = NULL;
}

SurfaceImpedanceBC::~SurfaceImpedanceBC() {

}

SurfaceImpedanceBC::SurfaceImpedanceBC(
        const CellTet<ORDER_N>* cell,
        const UInt face,
        const CellTet<ORDER_N>* cellD,
        const UInt faceD,
        const PMSurfaceSIBC* sibc) :
            BoundaryCondition(cell, face) {
    cellD_ = cellD;
    faceD_ = faceD;
    sibc_ = sibc;
}

SurfaceImpedanceBC& SurfaceImpedanceBC::operator=(
        const SurfaceImpedanceBC &rhs) {
    if (this == &rhs) {
        return *this;
    }
    BoundaryCondition::operator =(rhs);
    cellD_ = rhs.cellD_;
    faceD_ = rhs.faceD_;
    sibc_ = rhs.sibc_;
    return *this;
}

const CellTet<ORDER_N>* SurfaceImpedanceBC::getCellD() const {
    return cellD_;
}

UInt SurfaceImpedanceBC::getFaceD() const {
    return faceD_;
}

const PMPredefined* PhysicalModelBC::getCondition() const {
    return pm_;
}

const EMSourceBase* EMSourceBC::getCondition() const {
    return em_;
}

const PMSurfaceSIBC* SurfaceImpedanceBC::getCondition() const {
    return sibc_;
}

const CellTet<ORDER_N>* BoundaryCondition::getCell() const {
    return cell_;
}

UInt BoundaryCondition::getFace() const {
    return face_;
}

