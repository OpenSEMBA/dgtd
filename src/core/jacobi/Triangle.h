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
#ifndef CUDG3D_JACOBI_TRIANGLE_H_
#define CUDG3D_JACOBI_TRIANGLE_H_

#include "Line.h"
#include "math/matrix/Dynamic.h"
#include "math/util/SpaceGenerator.h"

namespace Cudg3d {
namespace Jacobi {

using namespace SEMBA;
using namespace Math;

template <size_t N>
class Triangle {
public:
    static const std::size_t faces = 3;
    static const std::size_t dimension = 2;
    static const std::size_t nfp = N + 1;
    static constexpr std::size_t np = (N+1)*(N+2)/2;

    Triangle();

    std::vector<CVecR2>   getGaussLobattoPoints()  const;
    DynMatR getVandermondeMatrix(    const std::vector<CVecR2>& rs) const;
    DynMatR getGradVandermondeMatrix(const std::vector<CVecR2>& rs) const;
    DynMatR getDifferentiationMatrix(const std::vector<CVecR2>& rs) const;
    DynMatR getLiftMatrix(           const std::vector<CVecR2>& rs) const;

    static std::vector<Real> evaluatePolynomialAt(
            const std::vector<CVecR2>& x,
            const std::pair<size_t,size_t> ij);

    static std::vector<Real> evaluateGradPolynomialAt(
            const std::vector<CVecR2>& x,
            const std::pair<size_t,size_t> ij);
private:
    static std::vector<Real> warpFactor_(const std::vector<Real>& rOut);
    static std::vector<CVecR2> rsToab_(const std::vector<CVecR2>& rs);
};

} /* namespace Jacobi */
} /* namespace DGTD */

#include "Triangle.hpp"

#endif
