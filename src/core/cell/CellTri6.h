//// OpenSEMBA
//// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
////                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
////                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
////                    Daniel Mateos Romero            (damarro@semba.guru)
////
//// This file is part of OpenSEMBA.
////
//// OpenSEMBA is free software: you can redistribute it and/or modify it under
//// the terms of the GNU Lesser General Public License as published by the Free
//// Software Foundation, either version 3 of the License, or (at your option)
//// any later version.
////
//// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
//// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//// details.
////
//// You should have received a copy of the GNU Lesser General Public License
//// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
//
//#ifndef CELLTRI6_H_
//#define CELLTRI6_H_
//
//#include "../../../common/geometry/Surface.h"
//#include "../../dgtd/core/CellTri.h"
//
//template <int TRI_N>
//class CellTri6 : public CellTri<TRI_N>, public Tri6 {
//	static const size_t np = (TRI_N+1) * (TRI_N+2) / 2;
//	static const size_t nfp = (TRI_N + 1);
//	static const size_t faces = 3;
//	static const size_t vertices = 3;
//public:
//	CellTri6(
//	 const Tri6& base);
//	virtual
//	~CellTri6();
//	CVecC3
//	 getRadiatedField(
//	  const CVecC3 electric[np],
//	  const CVecC3 magnetic[np],
//	  const double frequency,
//	  const pair<double, double> direction) const;
//};
//
//#include "../../dgtd/core/CellTri6.hpp"
//
//#endif /* CELLTRI6_H_ */
