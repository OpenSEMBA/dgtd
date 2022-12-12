#include "geometry/element/Tetrahedron4.h"
#include "math/function/Polynomial.h"
#include "Cell.h"

namespace SEMBA::dgtd::dg {

template <int N>
class Cell {
public:
    static const Math::Simplex::Tetrahedron<N> tet;
    static const std::size_t np = Math::Simplex::Tetrahedron<N>::np;

    using MatNpNp = Math::Matrix::Static<Math::Real, np, np>;

    const Geometry::Tet* base;
    size_t vmapP[faces][nfp]; // Node to Node of the contiguous element.
    Math::CVecR3 n[np];       // Lagrange's base functions nodes pos.

    Cell(const Geometry::Tet* base_, const PMGroup& pMGroup);
    virtual ~Cell() = default;
   
    size_t getNodeVertex(const size_t i) const;
    Geometry::ElemId getId() const { return base->getId(); }
    const Geometry::CoordR3* getV(size_t i) const { return base->getV(i); }
    array<MatNpNp, 3> getCMatrices() const;

    Math::CVecR3 getNode(const size_t i) const { return n[i]; }

    
    Cell*& getEtoEPointer(const int) const;

    Math::CVecR3 getSideNodePos(const size_t f, const size_t i) const;

    size_t getSideNode(size_t f, size_t i) const;

    virtual MatNpNp getMassMatrix() const;

    Math::Matrix::Static<double,np,np> getCMatrix(
        const size_t x,
        const Math::Matrix::Static<double, np, np>& invM,
        const Math::Matrix::Static<double, 4, 3> cJHat[Math::Simplex::Tetrahedron<1>::ncp],
        const Math::Matrix::Static<double, 4, 4> cJ[Math::Simplex::Tetrahedron<1>::ncp]) const;

    virtual Math::Matrix::Static<double,np,np> getMassMatrix(
        const double cJDet[Math::Simplex::Tetrahedron<1>::ncp]) const;

    virtual Math::Matrix::Static<double,np,np> getMassMatrixIntegratedWithScalar(
        const double cScalar[Math::Simplex::Tetrahedron<1>::ncp]) const;

    Math::Matrix::Static<double, np, nfp> getLIFTMatrix(
        const size_t s,
        const Math::Matrix::Static<double, np, nfp>& invM,
        const  double csd[Math::Simplex::Triangle<1>::ncp]) const;
private:
    void buildNodes();
};

template <int N>
const Math::Simplex::Tetrahedron<N> Cell<N>::tet;

template <int N>
Cell<N>::Cell(const Tet* base_, const PMGroup& pMGroup)
{
    base = base_;    
    material = pMGroup->getId(base->getMatId())->castTo<PMVolumeClassic>();
    buildNodes();
}

template <int N>
array<typename Cell<N>::MatNpNp, 3> Cell<N>::getCMatrices() const 
{
    MatR44 cJ[SimplexTet<1>::ncp];
    base->getCubatureJacobian(cJ);
    double cJDet[SimplexTet<1>::ncp];
    base->getCubatureJacobianDeterminant(cJDet, cJ);
    MatNpNp invM = getMassMatrix(cJDet).invert();
    StaMatrix<double,4,3> cJHat[SimplexTet<1>::ncp];
    base->getCubatureJacobianHat(cJHat, cJ, cJDet);
    array<MatNpNp,3> res;
    for (size_t x = 0; x < 3; x++) {
        res[x] = getCMatrix(x, invM, cJHat, cJ);
    }
    return res;
}

template <int N>
StaMatrix<double,np,np> Cell<N>::getCMatrix(
    const size_t x,
    const StaMatrix<double,np,np>& invM,
    const StaMatrix<double,4,3> cJHat[SimplexTet<1>::ncp],
    const StaMatrix<double,4,4> cJ[SimplexTet<1>::ncp]) const 
{
    StaMatrix<double,TET_NP,TET_NP> res;
    // Computes preliminary C matrices.
    for (size_t i = 0; i < this->tet.np; i++) {
        for (size_t j = 0; j < this->tet.np; j++) {
            for (size_t k = 0; k < faces; k++) {
                for (size_t c = 0; c < SimplexTet<1>::ncp; c++) {
                    res(i,j) += this->tet.cwada[c][k](i,j) * cJHat[c](k,x);
                }
            }
        }
    }
    // Multiplies by invMM.
    res = invM * res;
    res *= double(1.0 / 6.0);
    return res;
}

template <int N>
StaMatrix<double,TET_NP,TET_NP> Cell<N>::getMassMatrix(
        const double cJDet[SimplexTet<1>::ncp]) const
{
    StaMatrix<double,TET_NP,TET_NP> res;
    for (size_t c = 0; c < SimplexTet<1>::ncp; c++) {
        res += tet.cwaa[c] * cJDet[c];
    }
    res *= double(1.0 / 6.0);
    return res;
}

template <int N>
StaMatrix<double,TET_NP,TET_NP> Cell<N>::getMassMatrix() const 
{
    static const size_t ncp = SimplexTet<1>::ncp;
    MatR44 cJ[SimplexTet<1>::ncp];
    base->getCubatureJacobian(cJ);
    double cJDet[ncp];
    base->getCubatureJacobianDeterminant(cJDet, cJ);
    return getMassMatrix(cJDet);
}

template <int N>
StaMatrix<double,TET_NP,TET_NP>
Cell<N>::getMassMatrixIntegratedWithScalar(const double cScalar[SimplexTet<1>::ncp]) const 
{
    static const size_t ncp = SimplexTet<1>::ncp;
    MatR44 cJ[SimplexTet<1>::ncp];
    base->getCubatureJacobian(cJ);
    double cJDetByScalar[ncp];
    base->getCubatureJacobianDeterminant(cJDetByScalar, cJ);
    for (size_t c = 0; c < ncp; c++) {
        cJDetByScalar[c] *= cScalar[c];
    }
    return getMassMatrix(cJDetByScalar);
}

template <int N>
size_t Cell<N>::getNodeVertex(const size_t i) const 
{
    return tet.vertex(i);
}

template <int N>
CVecR3 Cell<N>::getSideNodePos(const size_t f, const size_t i) const 
{
    return n[tet.sideNode(f, i)];
}

template <int N>
void Cell<N>::buildNodes() 
{
    // Evaluates Lagrange's functions in positions specified by the
    // simplex coordinates of tet.
    double lagrEv[tet.np][base->numberOfCoordinates()];
    for (size_t j = 0; j < tet.np; j++) {
        for (size_t i = 0; i < base->numberOfCoordinates(); i++) {
            lagrEv[j][i]= base->getTet().getLagr(i).eval(tet.coordinate(j));
        }
    }
    // Computes nodes.
    for (size_t j = 0; j < tet.np; j++) {
        for (size_t i = 0; i < base->numberOfCoordinates(); i++) {
            this->n[j] += *(base->getV(i)) * lagrEv[j][i];
        }
    }
}

}