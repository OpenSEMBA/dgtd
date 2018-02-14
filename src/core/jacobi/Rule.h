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
#ifndef CUDG3D_JACOBI_RULE_H_
#define CUDG3D_JACOBI_RULE_H_

#include "math/util/Real.h"
#include "math/Constants.h"

#include <vector>
#include <utility>
#include <stdexcept>
#include <numeric>

namespace Cudg3d {
namespace Jacobi {

using namespace SEMBA::Math;

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
