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

#ifndef CELLTRI_H_
#define CELLTRI_H_

#include <complex>
#include <array>

#include "Cell.h"

namespace SEMBA {
namespace Cudg3d {
namespace Cell {

template<int TRI_N>
class Triangle : public Cell {
public:
    static const Math::Simplex::Triangle<TRI_N> tri;
    static const size_t np       = (TRI_N+1) * (TRI_N+2) / 2;
    static const size_t nfp      = (TRI_N + 1);

    Triangle();
    virtual ~Triangle();

//    virtual Math::CVecC3 getRadiatedField(
//            const Math::CVecC3 electric[np],
//            const Math::CVecC3 magnetic[np],
//            const double frequency,
//            const std::pair<double,double> direction) const = 0;

};

}
}
}


#include "Triangle.hpp"

#endif /* CELLTRI_H_ */
