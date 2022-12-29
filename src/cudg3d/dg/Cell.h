#include <Eigen/Dense>

#include "physicalModel/volume/Classic.h"
#include "geometry/element/Tetrahedron4.h"

namespace SEMBA::cudg3d::dg {

class Cell {
public:
    using StiffnessMat = Eigen::MatrixXd;
    using LiftMat = Eigen::MatrixXd;

    using PolynomialOrder = std::size_t;

    using Direction = std::size_t;

    using Node = Math::CVecR3;
    using NodeId = std::size_t;

    using Face = std::size_t;
    
    PolynomialOrder order{ 1 };
    const Geometry::Tet4& base;
    const PhysicalModel::Volume::Classic& material;

    Cell(const PolynomialOrder&, const Geometry::Tet4& base_, const PMGroup& pMGroup);
    virtual ~Cell() = default; 
    
    StiffnessMat getCMatrix(const Direction&) const;
    LiftMat getLiftMatrix(const Face&) const;

    std::size_t getNumberOfNodes() const;
    std::size_t getNumberOfFaceNodes() const;
};

}