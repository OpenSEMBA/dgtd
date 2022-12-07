#pragma once

#include "math/util/Real.h"
#include "math/Constants.h"

#include <vector>
#include <utility>
#include <stdexcept>
#include <numeric>

namespace SEMBA::dgtd::jacobi {

using namespace Math;

class Rule {
public:
    Rule(std::size_t order,
         std::pair<Real,Real> alphabeta,
         std::pair<Real,Real> region);

    const std::vector<Real>& getWeights() const {
        return w_;
    }

    const std::vector<Real>& getPoints() const {
        return x_;
    }

private:
    const size_t order_;
    const Real alpha_;
    const Real beta_;
    const Real a_;
    const Real b_;

    std::vector<Real> x_;
    std::vector<Real> w_;

    void cdgqf(std::vector<double>& t, std::vector<double>& wts) const;
    void cgqf(std::vector<double>& t, std::vector<double>& wts) const;
    double class_matrix (std::vector<double>& aj, std::vector<double>& bj) const;
    void imtqlx(
            std::vector<double>& d,
            std::vector<double>& e,
            std::vector<double>& z) const;
    void scqf(
            const std::vector<double>& t,
            const std::vector<size_t>& mlt,
            const std::vector<double>& wts,
            const std::vector<size_t>& ndx,
            std::vector<double>& swts,
            std::vector<double>& st) const ;
    void sgqf(
            const std::vector<double>& aj, std::vector<double>& bj,
            double zemu,
            std::vector<double>& t, std::vector<double>& wts) const ;

}; /* Class Rule */

} /* namespace Jacobi */
} /* namespace DGTD */

#endif
