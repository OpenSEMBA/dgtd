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

#include "Triangle.h"

namespace Cudg3d {
namespace Jacobi {

using namespace SEMBA;
using namespace Math;

template <size_t N>
Triangle<N>::Triangle() {
};

template <size_t N>
std::vector<CVecR2> Triangle<N>::getGaussLobattoPoints() const {
    std::array<CVecR3,np> L;
    size_t sk = 0;
    for (size_t n = 0; n < N; ++n) {
        for (size_t m = 0; m < N+1-n; ++m) {
            Real L1 = (Real) (n-1) / (Real) N;
            Real L3 = (Real) (m-1) / (Real) N;
            Real L2 = 1.0 - L1 - L3;
            L[sk++] = {L1, L2, L3};
        }
    }

    // Blending function.
    std::array<CVecR3,np> blend;
    for (size_t i = 0; i < blend.size(); ++i) {
        blend[i] = {4.0 * L[i][1] * L[i][2],
                    4.0 * L[i][2] * L[i][0],
                    4.0 * L[i][0] * L[i][1]};
    }

    // Warp factor.
    std::array<CVecR3,np> warpf;
    {
        std::array<std::vector<Real>,3> diff;
        for (size_t i = 0; i < L.size(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                diff[j].push_back(L[i][(j+2)%3] - L[i][(j+1)%3]);
            }
        }
        std::array<std::vector<Real>,3> warpfAux;
        for (size_t j = 0; j < 3; ++j) {
            warpfAux[j] = warpFactor_(diff[j]);
        }
        for (size_t i = 0; i < warpf.size(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                warpf[i][j] = warpfAux[j][i];
            }
        }
    }

    // Combine blend & warp
    static std::array<Real,15> alpopt = {
            0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
            0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
            1.4959, 1.5743, 1.5770, 1.6223, 1.6258 };
    Real alpha;
    N < 15? alpha = alpopt[N] : alpha = 5.0 / 3.0;

    std::array<CVecR3,np> warp;
    for (size_t i = 0; i < warp.size(); ++i) {
        for (size_t j = 0; j < 3; j++) {
            warp[i][j] = blend[i][j] *
                         warpf[i][j] *
                         std::pow(1.0 + alpha*L[i][j], 2);
        }
    }

    std::array<CVecR2,np> xy;
    for (size_t i = 0; i < np; i++) {
        xy[i] = {-L[i][1]+L[i][2],  (2.0*L[i][0]-L[i][1]-L[i][2]) / sqrt(3.0)};
        xy[i] += {warp[i][0] + std::cos(2.0/3.0*M_PI)*warp[i][1] +
                               std::cos(4.0/3.0*M_PI)*warp[i][2],
                               std::sin(2.0/3.0*M_PI)*warp[i][1] +
                               std::sin(4.0/3.0*M_PI)*warp[i][2] };
    }

    // Converts to coordinates in standard triangle.
    std::array<CVecR2,np> rs;
    for (size_t i = 0; i < np; i++) {
        Real L1, L2, L3;
        const Real &x = xy[i][0];
        const Real &y = xy[i][1];
        L1 = (std::sqrt(3.0)*y + 1.0) / 3.0;
        L2 = (-3.0*x - std::sqrt(3.0)*y + 2.0)/6.0;
        L3 = ( 3.0*x - std::sqrt(3.0)*y + 2.0)/6.0;
        rs[i][0] = -L2 + L3 - L1;
        rs[i][1] = -L2 - L3 + L1;
    }
    return rs;
}

template <size_t N>
std::vector<Real> Triangle<N>::warpFactor_(const std::vector<Real>& rOut) {
    Line<N> line;
    std::vector<Real> rEq  = Util::linspace({-1.0, 1.0}, N+1);
    Matrix::Dynamic<Real> Veq = line.getVandermondeMatrix(rEq);

    // Evaluate Lagrange polynomial at rOut
    DynMatR Pmat(N+1, rOut.size());
    for (size_t i = 0; i < N+1; ++i) {
        Pmat.cpToRow(i, line.evaluatePolynomialAt(rOut, 0, 0, i));
    }
    Veq.transpose();
    Pmat.invert();
    DynMatR Lmat = Veq * Pmat;
    Lmat.transpose();

    // Compute warp factor.
    DynMatR LGLrVec(line.getGaussLobattoPoints());
    DynMatR rEqVec(rEq);
    DynMatR warp = Lmat * (LGLrVec - rEqVec);

    // Scale factor.
    // TODO

    return warp.cpRowToVector(0);
}

template <size_t N>
std::vector<CVecR2> Triangle<N>::rsToab_(const std::vector<CVecR2>& rs) {
    std::vector<CVecR2> ab(rs.size());
    for (size_t i = 0; i < rs.size(); ++i) {
        if (rs[i](1) != 1.0) {
            ab[i](0) = 2.0 * (1.0 + rs[i](0)) / (1.0 - rs[i](1)) - 1.0;
        } else {
            ab[i](0) = -1.0;
        }
        ab[i](1) = rs[i](1);
    }
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getVandermondeMatrix(
        const std::vector<CVecR2>& rs) const {
    Matrix::Dynamic<Real> res(rs.size(), np);
    std::vector<CVecR2> ab = rsToab_(rs);
    size_t sk = 0;
    for (size_t i = 0; i < np; ++i) {
        for (size_t j = 0; j < np-i; ++j) {
            res.cpToCol(sk++, this->evaluatePolynomialAt(ab, {i,j}));
        }
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getGradVandermondeMatrix(
        const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> res(x.size(), N+1);
    for (size_t i = 0; i < (N+1); ++i) {
        res.cpToCol(i, this->evaluateGradPolynomialAt(x, alpha_, beta_, i));
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getDifferentiationMatrix(
        const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> Vr = this->getGradVandermondeMatrix(x);
    Matrix::Dynamic<Real> invV = this->getVandermondeMatrix(x).invert();
    return Vr * invV;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getLiftMatrix(
        const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> V = this->getVandermondeMatrix(x);
    Matrix::Dynamic<Real> Vtrans = this->getVandermondeMatrix(x).transpose();
    Matrix::Dynamic<Real> extractionMatrix(np,nfp*faces);
    extractionMatrix(0,0)     = 1.0;
    extractionMatrix(np-1, 1) = 1.0;
    return (V * Vtrans) * extractionMatrix;
}

template <size_t N>
std::vector<Real> Triangle<N>::evaluateGradPolynomialAt(
                const std::vector<Real>& x,
                const Real alpha,
                const Real beta,
                const size_t n) {
    if (n == 0) {
        return std::vector<Real>(x.size(), 0.0);;
    } else {
        auto res = evaluatePolynomialAt(x, alpha+1.0, beta+1.0, n-1);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] *= sqrt(n*(n + alpha + beta + 1));
        }
        return res;
    }
}

template <size_t N>
std::vector<Real> Triangle<N>::evaluatePolynomialAt(
        const std::vector<CVecR2>& ab,
        const std::pair<size_t,size_t> ij) {
    std::vector<Real> a(ab.size()), b(ab.size());
    for (size_t i = 0; i < ab.size(); i++) {
        a[i] = ab[i][0];
        b[i] = ab[i][1];
    }

    Line<N> line;
    std::vector<Real> h1 =
            line.evaluatePolynomialAt(a,          0.0, 0.0, ij.first );
    std::vector<Real> h2 =
            line.evaluatePolynomialAt(b, 2*ij.first+1, 0.0, ij.second);

    std::vector<Real> res(ab.size());
    for (size_t i = 0; i < ab.size(); ++i) {
        res[i] = sqrt(2.0) * h1[i] * h2[i] * pow(1.0 - b[i], i);
    }
    return res;
}

} /* namespace Jacobi */
} /* namespace DGTD */
