/*
 * CellTri.h
 *
 *  Created on: Jul 26, 2013
 *      Author: luis
 */

#ifndef CELLTRI_H_
#define CELLTRI_H_

#include "Cell.h"
#include <complex>
#include "math/Constants.h"
#include "math/MathUtils.h"
#include "math/SphericalVector.h"

template<int TRI_N>
class CellTri : public Cell {
public:
    static const SimplexTri<TRI_N> tri;
    static const size_t np = (TRI_N+1) * (TRI_N+2) / 2;
    static const size_t nfp = (TRI_N + 1);
    static const size_t faces = 3;
    static const size_t vertices = 3;
    CellTri();
    virtual
    ~CellTri();
    virtual CVecC3 getRadiatedField(
            const CVecC3 electric[np],
            const CVecC3 magnetic[np],
            const double frequency,
            const pair<double,double> direction) const = 0;
    virtual size_t getNumberOfVertices() const;
    virtual CVecR3 getSideNormal(const size_t s) const;
};

#include "../../dgtd/core/CellTri.hpp"

#endif /* CELLTRI_H_ */
