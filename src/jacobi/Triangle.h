#pragma once

#include "Line.h"
#include "math/matrix/Dynamic.h"
#include "math/util/SpaceGenerator.h"

#include <cmath>

namespace SEMBA::dgtd::jacobi {

using namespace Math;

using DynMatR = Matrix::Dynamic<Real>;

template <size_t N>
class Triangle {
public:
    static constexpr std::size_t faces = 3;
    static constexpr std::size_t dimension = 2;
    static constexpr std::size_t nfp = N + 1;
    static constexpr std::size_t np = (N+1)*(N+2)/2;
        
    std::vector<CVecR2>   getGaussLobattoPoints()  const;
    DynMatR getVandermondeMatrix(const std::vector<CVecR2>& rs) const;
    DynMatR getGradVandermondeMatrix(const std::vector<CVecR2>& rs) const;
    DynMatR getDifferentiationMatrix(const std::vector<CVecR2>& rs) const;
    DynMatR getLiftMatrix(const std::vector<CVecR2>& rs) const;

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

template <size_t N>
std::vector<CVecR2> Triangle<N>::getGaussLobattoPoints() const {
    std::array<CVecR3, np> L;
    size_t sk = 0;
    for (size_t n = 0; n < N; ++n) {
        for (size_t m = 0; m < N + 1 - n; ++m) {
            Real L1 = (Real)(n - 1) / (Real)N;
            Real L3 = (Real)(m - 1) / (Real)N;
            Real L2 = 1.0 - L1 - L3;
            L[sk++] = { L1, L2, L3 };
        }
    }

    // Blending function.
    std::array<CVecR3, np> blend;
    for (size_t i = 0; i < blend.size(); ++i) {
        blend[i] = { 4.0 * L[i][1] * L[i][2],
                    4.0 * L[i][2] * L[i][0],
                    4.0 * L[i][0] * L[i][1] };
    }

    // Warp factor.
    std::array<CVecR3, np> warpf;
    {
        std::array<std::vector<Real>, 3> diff;
        for (size_t i = 0; i < L.size(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                diff[j].push_back(L[i][(j + 2) % 3] - L[i][(j + 1) % 3]);
            }
        }
        std::array<std::vector<Real>, 3> warpfAux;
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
    static std::array<Real, 15> alpopt = {
            0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
            0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
            1.4959, 1.5743, 1.5770, 1.6223, 1.6258 };
    Real alpha;
    N < 15 ? alpha = alpopt[N] : alpha = 5.0 / 3.0;

    std::array<CVecR3, np> warp;
    for (size_t i = 0; i < warp.size(); ++i) {
        for (size_t j = 0; j < 3; j++) {
            warp[i][j] = blend[i][j] *
                warpf[i][j] *
                std::pow(1.0 + alpha * L[i][j], 2);
        }
    }

    std::array<CVecR2, np> xy;
    for (size_t i = 0; i < np; i++) {
        xy[i] = { -L[i][1] + L[i][2],  (2.0 * L[i][0] - L[i][1] - L[i][2]) / sqrt(3.0) };
        xy[i] += {warp[i][0] + std::cos(2.0 / 3.0 * M_PI) * warp[i][1] +
            std::cos(4.0 / 3.0 * M_PI) * warp[i][2],
            std::sin(2.0 / 3.0 * M_PI)* warp[i][1] +
            std::sin(4.0 / 3.0 * M_PI) * warp[i][2] };
    }

    // Converts to coordinates in standard triangle.
    std::vector<CVecR2> rs(np);
    for (size_t i = 0; i < np; i++) {
        Real L1, L2, L3;
        const Real& x = xy[i][0];
        const Real& y = xy[i][1];
        L1 = (std::sqrt(3.0) * y + 1.0) / 3.0;
        L2 = (-3.0 * x - std::sqrt(3.0) * y + 2.0) / 6.0;
        L3 = (3.0 * x - std::sqrt(3.0) * y + 2.0) / 6.0;
        rs[i][0] = -L2 + L3 - L1;
        rs[i][1] = -L2 - L3 + L1;
    }
    return rs;
}

template <size_t N>
std::vector<Real> Triangle<N>::warpFactor_(const std::vector<Real>& rOut) {
    Line<N> line;
    std::vector<Real> rEq = Util::linspace(
        std::pair<Real, Real>(-1.0, 1.0), N + 1);
    Matrix::Dynamic<Real> Veq = line.getVandermondeMatrix(rEq);

    // Evaluate Lagrange polynomial at rOut
    DynMatR Pmat(N + 1, rOut.size());
    for (size_t i = 0; i < N + 1; ++i) {
        Pmat.cpToRow(i, line.evaluatePolynomialAt(rOut, 0, 0, i));
    }
    Veq.transpose();
    Pmat.invert();
    DynMatR Lmat;
    Lmat = Veq * Pmat;
    Lmat.transpose();

    // Compute warp factor.
    DynMatR LGLrVec(rOut.size(), 1);
    LGLrVec.cpToCol(0, line.getGaussLobattoPoints());
    DynMatR rEqVec(rEq.size(), 1);
    rEqVec.cpToCol(0, rEq);
    DynMatR warp;
    warp = Lmat * (LGLrVec + rEqVec * (-1.0));

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
        }
        else {
            ab[i](0) = -1.0;
        }
        ab[i](1) = rs[i](1);
    }
    return ab;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getVandermondeMatrix(
    const std::vector<CVecR2>& rs) const {
    Matrix::Dynamic<Real> res(rs.size(), np);
    std::vector<CVecR2> ab = rsToab_(rs);
    size_t sk = 0;
    for (size_t i = 0; i < np; ++i) {
        for (size_t j = 0; j < np - i; ++j) {
            res.cpToCol(sk++, this->evaluatePolynomialAt(ab, { i,j }));
        }
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getGradVandermondeMatrix(
    const std::vector<CVecR2>& x) const {
    Matrix::Dynamic<Real> res(x.size(), N + 1);
    for (size_t i = 0; i < (N + 1); ++i) {
        res.cpToCol(i, this->evaluateGradPolynomialAt(x, i));
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getDifferentiationMatrix(
    const std::vector<CVecR2>& x) const {
    Matrix::Dynamic<Real> Vr = this->getGradVandermondeMatrix(x);
    Matrix::Dynamic<Real> invV = this->getVandermondeMatrix(x).invert();
    return Vr * invV;
}

template <size_t N>
Matrix::Dynamic<Real> Triangle<N>::getLiftMatrix(
    const std::vector<CVecR2>& x) const {
    Matrix::Dynamic<Real> V = this->getVandermondeMatrix(x);
    Matrix::Dynamic<Real> Vtrans = this->getVandermondeMatrix(x).transpose();
    Matrix::Dynamic<Real> extractionMatrix(np, nfp * faces);
    extractionMatrix(0, 0) = 1.0;
    extractionMatrix(np - 1, 1) = 1.0;
    return (V * Vtrans) * extractionMatrix;
}

template <size_t N>
std::vector<Real> Triangle<N>::evaluateGradPolynomialAt(
    const std::vector<CVecR2>& rs,
    const std::pair<size_t, size_t> ij) {
    //    if (n == 0) {
    //        return std::vector<Real>(rs.size(), 0.0);;
    //    } else {
    //        auto res = evaluatePolynomialAt(rs, n-1);
    //        for (size_t i = 0; i < res.size(); ++i) {
    //            res[i] *= sqrt(n*(n + 1));
    //        }
    //        return res;
    //    }
}

template <size_t N>
std::vector<Real> Triangle<N>::evaluatePolynomialAt(
    const std::vector<CVecR2>& ab,
    const std::pair<size_t, size_t> ij) {
    std::vector<Real> a(ab.size()), b(ab.size());
    for (size_t i = 0; i < ab.size(); i++) {
        a[i] = ab[i][0];
        b[i] = ab[i][1];
    }

    Line<N> line;
    auto h1 = line.evaluatePolynomialAt(a, 0.0, 0.0, ij.first);
    auto h2 = line.evaluatePolynomialAt(b, 2 * ij.first + 1, 0.0, ij.second);

    std::vector<Real> res(ab.size());
    for (size_t i = 0; i < ab.size(); ++i) {
        res[i] = sqrt(2.0) * h1[i] * h2[i] * pow(1.0 - b[i], i);
    }
    return res;
}


}