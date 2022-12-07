#pragma once 

#include "BoundaryCondition.h"
#include "mesh/Volume.h"
#include "source/Group.h"
#include "physicalModel/Group.h"
#include "geometry/graph/Connectivities.h"
#include "geometry/element/Hexahedron8.h"

namespace SEMBA {
namespace Cudg3d {
namespace BoundaryCondition {

template<typename T = Base>
class Group : public SEMBA::Group::Group<T> {
public:
    Group(const Mesh::Volume& mesh,
            const SourceGroup& em,
            const PMGroup& pm);
    void printInfo() const;

private:
    void buildEMSourceBC_(
            const Mesh::Volume& mesh,
            const SourceGroup& em);
    void buildPhysicalModelBC_(
            const Mesh::Volume& mesh,
            const PMGroup& pm);
    void removeOverlapped();
    vector<Base*> removeCommons(
            const vector<Base*>& low,
            const vector<Base*>& high) const;
};

template<typename T>
Group<T>::Group(
        const Mesh::Volume& mesh,
        const SourceGroup& em,
        const PMGroup& pm) {
    buildEMSourceBC_(mesh, em);
    buildPhysicalModelBC_(mesh, pm);
    removeOverlapped();
}

template<typename T>
void Group<T>::buildEMSourceBC_(
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
            add(new SourceBC(localFace, neighFace, em(i) ) );
        }
    }
}

template<typename T>
void Group<T>::buildPhysicalModelBC_(
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
            const PhysicalModel::PhysicalModel* mat = surf(i)->getModel();
            add(new PhysicalModelBC(localFace, neighFace, mat));
        }
    }
}