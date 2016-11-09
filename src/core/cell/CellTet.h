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

#ifndef CELLTET_H_
#define CELLTET_H_

#include <vector>
#include <set>

using namespace std;
using namespace SEMBA;

#include "math/simplex/Tetrahedron.h"
#include "math/function/Polynomial.h"
#include "Cell.h"

namespace SEMBA {
namespace Cudg3d {
namespace Cell {

template <int TET_N>
class CellTet : public Cell {
#define TET_NP ((TET_N+1)*(TET_N+2)*(TET_N+3)/6)
#define TET_NFP ((TET_N+1)*(TET_N+2)/2)
public:
    static const size_t np = TET_NP;
    static const size_t nfp = TET_NFP;
    static const size_t faces = 4;
    static const size_t vertices = 4;
    typedef Math::Matrix::Static<Math::Real,TET_NP,TET_NP> MatNpNp;
    static const Math::Simplex::Tetrahedron<TET_N> tet;
    const Geometry::Tet* base;
    size_t vmapP[faces][nfp]; // Node to Node of the contiguous element.
    Math::CVecR3 n[np]; // Lagrange's base functions nodes pos.

    CellTet();
    virtual ~CellTet();
    bool isCurved() const;
    bool isCurvedFace(const size_t f) const;
    double getVolume() const;
    double getAreaOfFace(size_t face) const;
    bool isLocalSide(
            const size_t side,
            const Geometry::SurfR* surf) const;
    bool isLocalSide(const Geometry::SurfR* surf) const;
    size_t getFaces() const {return faces;}
    size_t getNbp() const {return base->numberOfCoordinates();}
    size_t getNbfp() const;
    size_t getNfp() const {return nfp;}
    size_t getNumberOfVertices() const {return vertices;}
    size_t getNodeVertex(const size_t i) const;
    Geometry::ElemId getId() const {return base->getId();}
    const Geometry::CoordR3* getV(size_t i) const {return base->getV(i);}
    array<MatNpNp,3> getCMatrices() const;
    Math::CVecR3 getNode(const size_t i) const {return n[i];}
    const Geometry::CoordR3* getSideBaseNode(
            const size_t f,
            const size_t i) const;
    const Geometry::CoordR3* getSideVertexBaseNode(size_t f, size_t i) const;
    Cell*& getEtoEPointer(const int) const;
    Math::CVecR3 getSideNodePos(const size_t f,const size_t i) const;
    Math::CVecR3 getSideNormal(const size_t f) const;
    virtual void getCurvedLIFTnormal(
            Math::Matrix::Static<double,np,nfp> LIFTn[3],
            Math::Matrix::Static<double,np,nfp> LIFTcn[3],
            Math::Matrix::Static<double,np,nfp> LIFTrn[3],
            const size_t face) const;
    const Geometry::Tet* getBase() const;
    size_t getSideNode(size_t f, size_t i) const;
    bool isFaceContainedInPlane(
            const size_t face,
            const Math::Constants::CartesianPlane) const;
    virtual MatNpNp getConductivityWithGeometricProfile(
            const PhysicalModel::Volume::PML& mat,
            const size_t type,
            const double maxSigma) const;
    virtual MatNpNp getMassMatrix() const;
    virtual void printInfo() const;
    void printMapsInfo() const;
protected:
    void init(
            const Geometry::Tet* base_,
            const PMGroup* pMGroup);
    Math::Matrix::Static<double,TET_NP,TET_NP> getCMatrix(
            const size_t x,
            const Math::Matrix::Static<double,np,np>& invM,
            const Math::Matrix::Static<double,4,3> cJHat[Math::Simplex::Tetrahedron<1>::ncp],
            const Math::Matrix::Static<double,4,4> cJ[Math::Simplex::Tetrahedron<1>::ncp]) const;
    virtual Math::Matrix::Static<double,TET_NP,TET_NP> getMassMatrix(
            const double cJDet[Math::Simplex::Tetrahedron<1>::ncp]) const;
    virtual Math::Matrix::Static<double,TET_NP,TET_NP> getMassMatrixIntegratedWithScalar(
            const double cScalar[Math::Simplex::Tetrahedron<1>::ncp]) const;
    Math::Matrix::Static<double,np,nfp> getLIFTMatrix(
            const size_t s,
            const Math::Matrix::Static<double,np,nfp>& invM,
            const  double csd[Math::Simplex::Triangle<1>::ncp]) const;
private:
    void buildNodes();
};

template <int TET_N>
const Math::Simplex::Tetrahedron<TET_N> CellTet<TET_N>::tet;

template <int TET_N>
class CellTet4 : public CellTet<TET_N> {
    friend class Tetrahedron4;
public:
    static const size_t np = TET_NP;
    static const size_t nfp = TET_NFP;
    static const size_t faces = 4;
    static const size_t vertices = 4;
    CellTet4();
    virtual ~CellTet4();
    CellTet4(const Geometry::Tet* element, const PMGroup& pMGroup);
    bool isCurved() const;
    void printBCInfo() const;
    void printMapsInfo() const;
};

template <int TET_N>
class CellTet10 : public CellTet<TET_N> {
public:
    static const size_t np = TET_NP;
    static const size_t nfp = TET_NFP;
    static const size_t faces = 4;
    static const size_t vertices = 4;
    //
    CellTet10();
    virtual ~CellTet10();
    CellTet10(const Geometry::Tet* element, const PMGroup& pMGroup);
    void getCurvedLIFTnormal(
            Math::Matrix::Static<double,np,nfp> LIFTn[3],
            Math::Matrix::Static<double,np,nfp> LIFTcn[3],
            Math::Matrix::Static<double,np,nfp> LIFTrn[3],
            const size_t face) const;
    void buildOperators();
};

#include "CellTet.hpp"

}
}
}

#endif /* CELL_H_ */
