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
    points_.front() = -1.0;
    points_.back()  =  1.0;
    if (N != 1) {
        Rule rule(N-1,
                  std::make_pair(alpha_+1.0, beta_+1.0),
                  std::make_pair(-1.0, +1.0));

        const std::vector<Real> rulePoints = rule.getPoints();
        for (size_t i = 1; i < N; ++i) {
            points_[i] = rulePoints[i-1];
        }
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

} /* namespace Jacobi */
} /* namespace DGTD */
