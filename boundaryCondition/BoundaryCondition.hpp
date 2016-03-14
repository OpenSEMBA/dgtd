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

#include "BoundaryCondition.h"

namespace SEMBA {
namespace Cudg3d {
namespace BoundaryCondition {

template<class T>
BoundaryCondition<T>::BoundaryCondition::BoundaryCondition(
        T* condition,
        Cell::Face localFace,
        Cell::Face neighFace) {
    condition_ = condition;
    localFace_ = localFace;
    neighFace_ = neighFace;
}

template<class T>
BoundaryCondition<T>::~BoundaryCondition() {

}

template<class T>
bool BoundaryCondition<T>::hasSameBoundary(
        const BoundaryCondition<T>& other) const {
    return (getLocalFace() == other.getLocalFace());
}

template<class T>
BoundaryCondition<T>& BoundaryCondition<T>::operator =(
        const BoundaryCondition<T>& rhs) {
    if (this == &rhs) {
        return *this;
    }
    condition_ = rhs.condition_;
    localFace_ = rhs.localFace_;
    neighFace_ = rhs.neighFace_;
    return *this;
}

template<class T>
void BoundaryCondition<T>::printInfo() const {
    cout << "--- BC info ---" << endl;
    this->getCondition()->printInfo();

    cout << "Local face: " << endl;
    this->getLocalFace().first.printInfo();
    cout << "Face #: " << this->getLocalFace().second << endl;

    cout << "Neigh face: " << endl;
    this->getNeighFace().first.printInfo();
    cout << "Face #: " << this->getNeighFace().second << endl;
    cout << "--- End of BC Info ---" << endl;
}


template<class T>
inline Cell::Face BoundaryCondition<T>::getLocalFace() const {
    return localFace_;
}

template<class T>
inline Cell::Face BoundaryCondition<T>::getNeighFace() const {
    return neighFace_;
}

template<class T>
inline const T* BoundaryCondition<T>::getCondition() const {
    return condition_:
}

}
}
}

