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

#include "Line.h"

namespace Cudg3d {
namespace Jacobi {

using namespace SEMBA;
using namespace Math;

template <size_t N>
Line<N>::Line(Real alpha, Real beta) :
        alpha_(alpha),
        beta_(beta) {

    // Computes quadrature points.
    if (N == 1) {
        points_[0] = -1.0;
        points_[1] =  1.0;
    } else {
        jacobiGaussQuadrature(alpha_+1.0, beta_+1.0);
    }
};

template <size_t N>
inline std::vector<typename Line<N>::Point> Line<N>::getPoints() const {
    std::vector<Point> res(np);
    std::copy_n(points_.begin(), np, res.begin());
    return res;
}

template <size_t N>
inline std::vector<Real> Line<N>::getWeights() const {
    std::vector<Weight> res(np);
    std::copy_n(weights_.begin(), np, res.begin());
    return res;
}

template <size_t N>
std::vector<std::pair<typename Line<N>::Point, typename Line<N>::Weight>>
        Line<N>::jacobiGaussQuadrature(Real alpha, Real beta) {

    size_t constexpr n = N -2;
    std::static_assert(n >= 0);

    std::vector<std::pair<Point, Weight>> res;
    if (n == 0) {
        Point x = (alpha - beta)/(alpha + beta + 2.0);
        Weight w = 2.0;
        res.push_back(std::make_pair(x,w));
        return res;
    }

    Matrix::Static<Real, n+1, n+1> J;

#error "TO DO." // TODO
}

} /* namespace Jacobi */
} /* namespace DGTD */
