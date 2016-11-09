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
#ifndef CELL_TRIANGLE3_H_
#define CELL_TRIANGLE3_H_

#include "geometry/element/Triangle3.h"
#include "Triangle.h"

namespace SEMBA {
namespace Cudg3d {
namespace Cell {

template<int N>
class Triangle3 : public Triangle<N>, public Geometry::Tri3 {
public:
	Triangle3(const Geometry::Tri3& base);
	virtual	~Triangle3();

//	Math::CVecC3 getRadiatedField(
//	        const CVecC3 electric[np],
//	        const CVecC3 magnetic[np],
//	        const double frequency,
//	        const std::pair<double, double> direction) const;
};

}  /* namespace Cell */
}  /* namespace Cudg3d */
}  /* namespace SEMBA */

#include "Triangle3.hpp"

#endif /* CELL_TRIANGLE3_H_ */
