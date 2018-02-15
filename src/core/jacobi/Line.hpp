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
    gaussLobattoPoints_.front() = -1.0;
    gaussLobattoPoints_.back()  =  1.0;
    if (N != 1) {
        Rule rule(N-1,
                  std::make_pair(alpha_+1.0, beta_+1.0),
                  std::make_pair(-1.0, +1.0));

        const std::vector<Real> rulePoints = rule.getPoints();
        for (size_t i = 1; i < N; ++i) {
            gaussLobattoPoints_[i] = rulePoints[i-1];
        }
    }
};

template <size_t N>
inline std::vector<Real> Line<N>::getGaussLobattoPoints() const {
    std::vector<Real> res(np);
    std::copy_n(gaussLobattoPoints_.begin(), np, res.begin());
    return res;
}

template <size_t N>
std::vector<Real>  Line<N>::evaluateAt(const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> PL(np, x.size());
    const Real gamma0 = std::pow(2.0, alpha_+beta_+1.0) /
            (alpha_ +beta_ + 1.0) *
            std::tgamma(alpha_ + 1.0) *
            std::tgamma(beta_  + 1.0) /
            std::tgamma(alpha_+beta_+1.0);
    for (size_t i = 0; i < PL.nCols(); ++i) {
        PL(0,i) = 1.0 / sqrt(gamma0);
    }
    if (N != 0) {
        const Real gamma1 = (alpha_+1.0) * (beta_+1.0) /
                            (alpha_ + beta_ + 3.0) *
                            gamma0;
        for (size_t i = 0; i < PL.nCols(); ++i) {
            PL(1,i) =  ( (alpha_+beta_+2.0)*x[i] / 2.0+(alpha_-beta_)/2.0 ) /
                        sqrt(gamma1);
        }


        Real aOld = 2.0 /
                    (2.0 + alpha_ + beta_) *
                    sqrt((alpha_+1.0)*(beta_+1.0) / (alpha_+beta_+3.0));
        for (size_t i = 0; i < (N-1); ++i) {
            const Real h1 = (Real) (2*(i+1)) + alpha_ + beta_;
            const Real aNew = 2.0 / (h1 + 2.0) *
                             sqrt( ( (Real) i + 2.0 ) *
                                   ( (Real) i + 2.0 + alpha_ + beta_)*
                                   ( (Real) i + 2.0 + alpha_) *
                                   ( (Real) i + 2.0 + beta_) /
                                   (h1 + 1.0) /
                                   (h1 + 3.0)
                                 );
             const Real bNew = - (pow(alpha_,2) - pow(beta_,2)) /
                                  h1 /
                                  (h1 + 2.0);
             for (size_t j = 0; j < PL.nCols(); ++j) {
                 PL(i+2,j) = 1.0 / aNew *
                             (-aOld*PL(i,j) + (x[j]-bNew) * PL(i+1,j) );
             }
            aOld = aNew;
        }
    }

    return PL.cpRowToVector(N);
}

} /* namespace Jacobi */
} /* namespace DGTD */
