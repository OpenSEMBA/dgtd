#include "Rule.h"

namespace SEMBA::dgtd::jacobi {

using namespace Math;

//   Copied from Burkardt:
//   http://people.sc.fsu.edu/~jburkardt/cpp_src/jacobi_rule/jacobi_rule.cpp
Rule::Rule(std::size_t order,
        std::pair<Real,Real> alphabeta,
        std::pair<Real,Real> region) :
                    order_(order),
                    alpha_(alphabeta.first),
                    beta_(alphabeta.second),
                    a_(region.first),
                    b_(region.second) {

    if (order < 0) {
        throw std::logic_error("Jacobi::Rule: Order must be >= 1.\n");
    }
    if (beta_ <= -1.0 ) {
        throw std::logic_error("Jacobi::Rule: BETA <= -1.0.\n");
    }

    x_.resize(order);
    w_.resize(order);
    cgqf(x_, w_);
}

void Rule::cdgqf(std::vector<double>& t, std::vector<double>& wts) const {
    std::vector<double> aj(order_), bj(order_);
    double zemu = class_matrix (aj, bj);
    sgqf(aj, bj, zemu, t, wts);
}

void Rule::cgqf (std::vector<double>& t, std::vector<double>& wts) const {
    cdgqf (t, wts );

    std::vector<size_t> mlt(order_, 1);
    std::vector<size_t> ndx(order_);
    std::iota(ndx.begin(), ndx.end(), 1);

    scqf (t, mlt, wts, ndx, wts, t);
}

double Rule::class_matrix(
        std::vector<double>& aj, std::vector<double>& bj) const {
    const double ab = alpha_ + beta_;
    const double a2b2 = beta_ * beta_ - alpha_ * alpha_;

    double abi = 2.0 + ab;
    aj[0] = (beta_ - alpha_) / abi;
    bj[0] = sqrt ( 4.0 * ( 1.0 + alpha_ ) * ( 1.0 + beta_ )
            / ( ( abi + 1.0 ) * abi * abi ) );
    for (size_t i = 2; i <= order_; i++ ) {
        abi = 2.0 * i + ab;
        aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
        abi = abi * abi;
        bj[i-1] = sqrt ( 4.0 * i * ( i + alpha_ ) * ( i + beta_ ) * ( i + ab )
                / ( ( abi - 1.0 ) * abi ) );
    }
    return pow(2.0, ab+1.0) *
            tgamma(alpha_+1.0) * tgamma(beta_+1.0) / tgamma (2.0 + ab);
}

void Rule::imtqlx(
        std::vector<double>& d,
        std::vector<double>& e,
        std::vector<double>& z) const {

    const size_t numberOfIterationsLimit = 30;

    if (order_ == 1) {
        return;
    }

    e[order_-1] = 0.0;
    for (size_t l = 1; l <= order_; l++ ) {
        size_t j = 0;
        while (true) {
            size_t m;
            for (m = l; m <= order_; m++ ) {
                if ( m == order_ ) {
                    break;
                }
                if (fabs(e[m-1]) <= std::numeric_limits<double>::epsilon()
                                    * (fabs(d[m-1]) + fabs(d[m]))) {
                    break;
                }
            }
            double p = d[l-1];
            if ( m == l ) {
                break;
            }
            if (j >= numberOfIterationsLimit) {
                throw std::logic_error(
                        "Jacobi::Rule: Iteration limit exceeded\n");
            }
            j = j + 1;
            double g = ( d[l] - p ) / ( 2.0 * e[l-1] );
            const double r =  sqrt ( g * g + 1.0 );
            g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * Util::sign(g) );
            double s = 1.0;
            double c = 1.0;
            p = 0.0;
            for (size_t ii = 1; ii <= (m-l); ii++ ) {
                size_t i = m - ii;
                double f = s * e[i-1];
                const double b = c * e[i-1];
                if ( fabs ( g ) <= fabs ( f ) ) {
                    c = g / f;
                    const double r =  sqrt ( c * c + 1.0 );
                    e[i] = f * r;
                    s = 1.0 / r;
                    c = c * s;
                } else {
                    s = f / g;
                    const double r =  sqrt ( s * s + 1.0 );
                    e[i] = g * r;
                    c = 1.0 / r;
                    s = s * c;
                }
                g = d[i] - p;
                const double r = ( d[i-1] - g ) * s + 2.0 * c * b;
                p = s * r;
                d[i] = g + p;
                g = c * r - b;
                f = z[i];
                z[i] = s * z[i-1] + c * f;
                z[i-1] = c * z[i-1] - s * f;
            }
            d[l-1] = d[l-1] - p;
            e[l-1] = g;
            e[m-1] = 0.0;
        }
    }

    // Sorts positions preseving order of associated weight.
    for (size_t ii = 2; ii <= order_; ii++ ) {
        const size_t i = ii - 1;
        size_t k = i;
        double p = d[i-1];
        for (size_t j = ii; j <= order_; j++ ) {
            if ( d[j-1] < p ) {
                k = j;
                p = d[j-1];
            }
        }
        if ( k != i ) {
            d[k-1] = d[i-1];
            d[i-1] = p;
            p = z[i-1];
            z[i-1] = z[k-1];
            z[k-1] = p;
        }
    }
}

void Rule::scqf (
        const std::vector<double>& t,
        const std::vector<size_t>& mlt,
        const std::vector<double>& wts,
        const std::vector<size_t>& ndx,
        std::vector<double>& swts,
        std::vector<double>& st) const {
    if (fabs(b_-a_) <= std::numeric_limits<double>::epsilon()) {
        throw std::logic_error("Jacobi::Rule: B - A too small.\n");
    }
    const double shft = ( a_ + b_ ) / 2.0;
    const double slp  = ( b_ - a_ ) / 2.0;
    const double p = pow ( slp, alpha_ + beta_ + 1.0 );
    for (size_t k = 0; k < order_; k++ ) {
        st[k] = shft + slp * t[k];
        size_t l = ndx[k];
        if ( l != 0 ) {
            double tmp = p;
            for (size_t i = l - 1; i <= l - 1 + mlt[k] - 1; i++ ) {
                swts[i] = wts[i] * tmp;
                tmp *= slp;
            }
        }
    }
}

void Rule::sgqf(
        const std::vector<double>& aj, std::vector<double>& bj,
        double zemu,
        std::vector<double>& t, std::vector<double>& wts) const {
    if (zemu <= 0.0){
        throw std::logic_error("Jacobi::Rule: ZEMU <= 0.\n");
    }
    t = aj;
    wts[0] = sqrt(zemu);
    for (size_t i = 1; i < order_; i++ ) {
        wts[i] = 0.0;
    }
    imtqlx (t, bj, wts);
    for (size_t i = 0; i < order_; i++ ) {
        wts[i] *= wts[i];
    }
}

}