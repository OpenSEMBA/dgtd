#include "Cell.h"

namespace SEMBA::cudg3d::dg {

Cell::Cell(const PolynomialOrder& n, const Geometry::Tet4& b, const PMGroup& pMGroup) :
    order{n},
    base{b},
    material{ *pMGroup.getId(base.getMatId())->castTo<PhysicalModel::Volume::Classic>() }
{
    if (!pMGroup.getId(base.getMatId())->is<PhysicalModel::Volume::Classic>()) {
        throw std::runtime_error("Cell must be made of a volumic material.");
    }
}

Cell::StiffnessMat Cell::getCMatrix(const Direction& x) const 
{
    StiffnessMat res;
    return res;
}

Cell::LiftMat Cell::getLiftMatrix(const Face& x) const
{
    LiftMat res;
    return res;
}

std::size_t Cell::getNumberOfNodes() const
{
    return 0;
}

std::size_t Cell::getNumberOfFaceNodes() const
{
    return 0;
}

}