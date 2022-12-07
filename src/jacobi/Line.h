#pragma once

#include <array>
#include <algorithm>
#include <utility>
#include <type_traits>

#include "math/matrix/Static.h"
#include "math/function/Polynomial.h"
#include "Rule.h"

namespace SEMBA::dgtd::jacobi {

using namespace Math;

template <size_t N>
class Line {
public:
    static const std::size_t faces = 2;
    static const std::size_t dimension = 1;
    static const std::size_t nfp = 1;
    static constexpr std::size_t np = N + 1;

    Line(Real alpha = 0.0, Real beta = 0.0);

    std::vector<Real> getGaussLobattoPoints()  const;
    Matrix::Dynamic<Real> getVandermondeMatrix(const std::vector<Real>& x) const;
    Matrix::Dynamic<Real> getGradVandermondeMatrix(const std::vector<Real>& x) const;
    Matrix::Dynamic<Real> getDifferentiationMatrix(
            const std::vector<Real>& x) const;
    Matrix::Dynamic<Real> getLiftMatrix(const std::vector<Real>& x) const;

    static std::vector<Real> evaluatePolynomialAt(
            const std::vector<Real>& x,
            const Real alpha,
            const Real beta,
            const size_t n);

    static std::vector<Real> evaluateGradPolynomialAt(
                const std::vector<Real>& x,
                const Real alpha,
                const Real beta,
                const size_t n);
private:
    Real alpha_, beta_;
};

template <size_t N>
Line<N>::Line(Real alpha, Real beta) :
    alpha_(alpha),
    beta_(beta) {
};

template <size_t N>
std::vector<Real> Line<N>::getGaussLobattoPoints() const {
    std::vector<Real> res(np);
    res.front() = -1.0;
    res.back() = 1.0;
    if (N != 1) {
        Rule rule(N - 1,
            std::make_pair(alpha_ + 1.0, beta_ + 1.0),
            std::make_pair(-1.0, +1.0));

        const std::vector<Real> rulePoints = rule.getPoints();
        for (size_t i = 1; i < N; ++i) {
            res[i] = rulePoints[i - 1];
        }
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Line<N>::getVandermondeMatrix(
    const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> res(x.size(), N + 1);
    for (size_t i = 0; i < (N + 1); ++i) {
        res.cpToCol(i, this->evaluatePolynomialAt(x, alpha_, beta_, i));
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Line<N>::getGradVandermondeMatrix(
    const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> res(x.size(), np);
    for (size_t i = 0; i < np; ++i) {
        res.cpToCol(i, this->evaluateGradPolynomialAt(x, alpha_, beta_, i));
    }
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Line<N>::getDifferentiationMatrix(
    const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> Vr = this->getGradVandermondeMatrix(x);
    Matrix::Dynamic<Real> invV = this->getVandermondeMatrix(x).invert();
    Matrix::Dynamic<Real> res;
    res = Vr * invV;
    return res;
}

template <size_t N>
Matrix::Dynamic<Real> Line<N>::getLiftMatrix(
    const std::vector<Real>& x) const {
    Matrix::Dynamic<Real> V = this->getVandermondeMatrix(x);
    Matrix::Dynamic<Real> Vtrans = this->getVandermondeMatrix(x).transpose();
    Matrix::Dynamic<Real> extractionMatrix(np, nfp * faces);
    extractionMatrix(0, 0) = 1.0;
    extractionMatrix(np - 1, 1) = 1.0;
    return (V * Vtrans) * extractionMatrix;
}

template <size_t N>
std::vector<Real> Line<N>::evaluateGradPolynomialAt(
    const std::vector<Real>& x,
    const Real alpha,
    const Real beta,
    const size_t n) {
    if (n == 0) {
        return std::vector<Real>(x.size(), 0.0);;
    }
    else {
        auto res = evaluatePolynomialAt(x, alpha + 1.0, beta + 1.0, n - 1);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] *= sqrt(n * (n + alpha + beta + 1));
        }
        return res;
    }
}

template <size_t N>
std::vector<Real> Line<N>::evaluatePolynomialAt(
    const std::vector<Real>& x,
    const Real alpha,
    const Real beta,
    const size_t n) {
    Matrix::Dynamic<Real> PL(n + 1, x.size());
    const Real gamma0 = std::pow(2.0, alpha + beta + 1.0) /
        (alpha + beta + 1.0) *
        std::tgamma(alpha + 1.0) *
        std::tgamma(beta + 1.0) /
        std::tgamma(alpha + beta + 1.0);
    for (size_t i = 0; i < PL.nCols(); ++i) {
        PL(0, i) = 1.0 / sqrt(gamma0);
    }
    if (n != 0) {
        const Real gamma1 = (alpha + 1.0) * (beta + 1.0) /
            (alpha + beta + 3.0) *
            gamma0;
        for (size_t i = 0; i < PL.nCols(); ++i) {
            PL(1, i) = ((alpha + beta + 2.0) * x[i] / 2.0 + (alpha - beta) / 2.0) /
                sqrt(gamma1);
        }


        Real aOld = 2.0 /
            (2.0 + alpha + beta) *
            sqrt((alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0));
        for (size_t i = 0; i < (n - 1); ++i) {
            const Real h1 = (Real)(2 * (i + 1)) + alpha + beta;
            const Real aNew = 2.0 / (h1 + 2.0) *
                sqrt(((Real)i + 2.0) *
                    ((Real)i + 2.0 + alpha + beta) *
                    ((Real)i + 2.0 + alpha) *
                    ((Real)i + 2.0 + beta) /
                    (h1 + 1.0) /
                    (h1 + 3.0)
                );
            const Real bNew = -(pow(alpha, 2) - pow(beta, 2)) /
                h1 /
                (h1 + 2.0);
            for (size_t j = 0; j < PL.nCols(); ++j) {
                PL(i + 2, j) = 1.0 / aNew *
                    (-aOld * PL(i, j) + (x[j] - bNew) * PL(i + 1, j));
            }
            aOld = aNew;
        }
    }

    return PL.cpRowToVector(n);
}

}