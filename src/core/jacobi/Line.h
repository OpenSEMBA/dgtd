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
#ifndef CUDG3D_JACOBI_LINE_H_
#define CUDG3D_JACOBI_LINE_H_

#include "math/matrix/Static.h"
#include "math/function/Polynomial.h"

#include <array>
#include <algorithm>
#include <pair>
#include <type_traits>

namespace Cudg3d {
namespace Jacobi {

using namespace SEMBA;
using namespace Math;

template <size_t N>
class Line {
public:
    static const std::size_t faces = 2;
    static const std::size_t dimension = 1;
    static const std::size_t nfp = 1;
    static constexpr std::size_t np = N + 1;

    typedef Real Point;
    typedef Real Weight;

    Line(Real alpha = 0.0, Real beta = 0.0);

    std::vector<Point>  getPoints()  const;
    std::vector<Weight> getWeights() const;

private:
    Real alpha_, beta_;
    std::array<Point,np>  points_;
    std::array<Weight,np> weights_;

    static std::vector<std::pair<Point,Weight>>
        jacobiGaussQuadrature(Real alpha, Real beta);
};

} /* namespace Jacobi */
} /* namespace DGTD */

#include "Line.hpp"

#endif
