#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>			// Stream I/O.
#include <cmath>
#include <vector>
#include <utility>

#include "geometry/element/Element.h"
#include "physicalModel/PhysicalModel.h"
#include "source/Source.h"

namespace SEMBA {
namespace dgtd {
namespace BoundaryCondition {
//
//class Base {
//public:
//    Base();
//    virtual ~Base();
//
//    virtual bool hasSameBoundary(const Base& other) const = 0;
//    virtual Geometry::Element::Face getLocalFace() const = 0;
//    virtual Geometry::Element::Face getNeighFace() const = 0;
//};
//
//template<class T>
//class BoundaryCondition : public Base,
//                          public virtual Class::Class,
//                          public virtual Class::Cloneable,
//                          public virtual Class::Shareable,
//                          public virtual Class::Printable{
//public:
//    BoundaryCondition(
//            T* condition,
//            Geometry::Element::Face localFace,
//            Geometry::Element::Face neighFace);
//    virtual ~BoundaryCondition();
//    bool hasSameBoundary(const BoundaryCondition::Base& other) const;
//    virtual BoundaryCondition& operator=(const BoundaryCondition& rhs);
//
//    Geometry::Element::Face getLocalFace() const;
//    Geometry::Element::Face getNeighFace() const;
//    const T* getCondition() const;
//
//    void printInfo() const;
//private:
//    const T* condition_;
//    Geometry::Element::Face localFace_, neighFace_;
//};
//
//
//BoundaryCondition::BoundaryCondition(
//        T* condition,
//        Geometry::Element::Face localFace,
//        Geometry::Element::Face neighFace) {
//    condition_ = condition;
//    localFace_ = localFace;
//    neighFace_ = neighFace;
//}
//
//BoundaryCondition::~BoundaryCondition() {
//
//}
//
//bool BoundaryCondition::hasSameBoundary(const Base& other) const {
//    return (getLocalFace() == other.getLocalFace());
//}
//
//BoundaryCondition& BoundaryCondition::operator =(
//        const BoundaryCondition& rhs) {
//    if (this == &rhs) {
//        return *this;
//    }
//    condition_ = rhs.condition_;
//    localFace_ = rhs.localFace_;
//    neighFace_ = rhs.neighFace_;
//    return *this;
//}
//
//template<class T>
//void BoundaryCondition<T>::printInfo() const {
//    cout << "--- BC info ---" << endl;
//    this->getCondition()->printInfo();
//
//    cout << "Local face: " << endl;
//    this->getLocalFace().first.printInfo();
//    cout << "Face #: " << this->getLocalFace().second << endl;
//
//    cout << "Neigh face: " << endl;
//    this->getNeighFace().first.printInfo();
//    cout << "Face #: " << this->getNeighFace().second << endl;
//    cout << "--- End of BC Info ---" << endl;
//}
//
//
//template<class T>
//inline Geometry::Element::Face BoundaryCondition<T>::getLocalFace() const {
//    return localFace_;
//}
//
//template<class T>
//inline Geometry::Element::Face BoundaryCondition<T>::getNeighFace() const {
//    return neighFace_;
//}
//
//template<class T>
//inline const T* BoundaryCondition<T>::getCondition() const {
//    return condition_;
//}

//typedef BoundaryCondition<Source::Base> SourceBC;
//typedef BoundaryCondition<PhysicalModel::PhysicalModel> PhysicalModelBC;

}
}
}

