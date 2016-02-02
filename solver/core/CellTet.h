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

#include "math/Simplex.h"
#include "math/MathMatrix.h"
#include "math/FunctionPolynomial.h"
#include "SmbData.h"
#include "Cell.h"

template <int TET_N>
class CellTet : public Cell {
#define TET_NP ((TET_N+1)*(TET_N+2)*(TET_N+3)/6)
#define TET_NFP ((TET_N+1)*(TET_N+2)/2)
public:
    typedef StaMatrix<Real,TET_NP,TET_NP> MatNpNp;
    static const SimplexTet<TET_N> tet;
    static const UInt np = TET_NP;
    static const UInt nfp = TET_NFP;
    static const UInt faces = 4;
    static const UInt vertices = 4;
    const Tet* base;
    UInt vmapP[faces][nfp]; // Node to Node of the contiguous element.
    CVecR3 n[np]; // Lagrange's base functions nodes pos.
    CellTet();
    virtual ~CellTet();
    bool isCurved() const;
    bool isCurvedFace(const UInt f) const;
    double getVolume() const;
    double getAreaOfFace(UInt face) const;
    bool isLocalSide(
            const UInt side,
            const SurfR* surf) const;
    bool isLocalSide(
            const SurfR* surf) const;
    UInt getFaces() const {return faces;}
    UInt getNbp() const {return base->numberOfCoordinates();}
    UInt getNbfp() const;
    UInt getNfp() const {return nfp;}
    UInt getNumberOfVertices() const {return vertices;}
    UInt getNodeVertex(const UInt i) const;
    ElementId getId() const {return base->getId();}
    const CoordR3* getV(UInt i) const {return base->getV(i);}
    array<MatNpNp,3> getCMatrices() const;
    CVecR3 getNode(const UInt i) const {return n[i];}
    const CoordR3* getSideBaseNode(
            const UInt f,
            const UInt i) const;
    const CoordR3* getSideVertexBaseNode(UInt f, UInt i) const;
    Cell*& getEtoEPointer(const int) const;
    CVecR3 getSideNodePos(const UInt f,const UInt i) const;
    CVecR3 getSideNormal(const UInt f) const;
    virtual void getCurvedLIFTnormal(
            StaMatrix<double,np,nfp> LIFTn[3],
            StaMatrix<double,np,nfp> LIFTcn[3],
            StaMatrix<double,np,nfp> LIFTrn[3],
            const UInt face) const;
    const Tet* getBase() const;
    UInt getSideNode(UInt f, UInt i) const;
    bool isFaceContainedInPlane(
            const UInt face,
            const CartesianPlane) const;
    virtual MatNpNp getConductivityWithGeometricProfile(
            const PMVolumePML& mat,
            const UInt type,
            const double maxSigma) const;
    virtual MatNpNp getMassMatrix() const;
    virtual void printInfo() const;
    void printMapsInfo() const;
protected:
    void init(
            const Tet* base_,
            const PMGroup* pMGroup);
    StaMatrix<double,TET_NP,TET_NP> getCMatrix(
            const UInt x,
            const StaMatrix<double,np,np>& invM,
            const StaMatrix<double,4,3> cJHat[SimplexTet<1>::ncp],
            const StaMatrix<double,4,4> cJ[SimplexTet<1>::ncp]) const;
    virtual StaMatrix<double,TET_NP,TET_NP> getMassMatrix(
            const double cJDet[SimplexTet<1>::ncp]) const;
    virtual StaMatrix<double,TET_NP,TET_NP> getMassMatrixIntegratedWithScalar(
            const double cScalar[SimplexTet<1>::ncp]) const;
    StaMatrix<double,np,nfp> getLIFTMatrix(
            const UInt s,
            const StaMatrix<double,np,nfp>& invM,
            const  double csd[SimplexTri<1>::ncp]) const;
private:
    void buildNodes();
};

template <int TET_N>
const SimplexTet<TET_N> CellTet<TET_N>::tet;

template <int TET_N>
class CellTet4 : public CellTet<TET_N> {
    friend class Tetrahedron4;
public:
    static const UInt np = TET_NP;
    static const UInt nfp = TET_NFP;
    static const UInt faces = 4;
    static const UInt vertices = 4;
    CellTet4();
    virtual ~CellTet4();
    CellTet4(const Tet* element, const PMGroup& pMGroup);
    bool isCurved() const;
    void printBCInfo() const;
    void printMapsInfo() const;
};

template <int TET_N>
class CellTet10 : public CellTet<TET_N> {
public:
    static const UInt np = TET_NP;
    static const UInt nfp = TET_NFP;
    static const UInt faces = 4;
    static const UInt vertices = 4;
    //
    CellTet10();
    virtual ~CellTet10();
    CellTet10(const Tet* element, const PMGroup& pMGroup);
    void getCurvedLIFTnormal(
            StaMatrix<double,np,nfp> LIFTn[3],
            StaMatrix<double,np,nfp> LIFTcn[3],
            StaMatrix<double,np,nfp> LIFTrn[3],
            const UInt face) const;
    void buildOperators();
};

#include "CellTet.hpp"

#endif /* CELL_H_ */
