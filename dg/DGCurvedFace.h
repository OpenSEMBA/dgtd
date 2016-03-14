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
/*
 * SolverCurvedFace.h
 *
 *  Created on: Mar 21, 2013
 *      Author: luis
 */

#ifndef DGCURVEDFACE_H_
#define DGCURVEDFACE_H_

#include "solver/dgtd/core/CellTet.h"
#include "math/MathMatrix.h"
#include "math/Field.h"

// TODO: Curved faces are buggy. An static remanent is left.
// TODO: Observed instability for extreme curvature.

class DGCurvedFace {
public:
	size_t solverPosition;
	DGCurvedFace();
	DGCurvedFace(
	 const Cell* cell,
	 const size_t f,
	 const size_t solverPosition,
	 FieldR3& rhsE_, FieldR3& rhsH_,
	 const FieldR3& dE_, const FieldR3& dH_,
	 const FieldR3& dresE_, const FieldR3& dresH_,
	 const Real impP,
	 const Real admP,
	 const Real impAv,
	 const Real admAv);
	virtual ~DGCurvedFace();
	DGCurvedFace&
	 operator=(const DGCurvedFace& rhs);
	void
	 addFluxToRHSElectric(
	  const Real upwinding,
	  const bool useResForUpw);
	void
	 addFluxToRHSMagnetic(
	  const Real upwinding,
	  const bool useResForUpw);
private:
	static const size_t N = ORDER_N;
	static const size_t np = (N+1) * (N+2) * (N+3) / 6;
	static const size_t nfp = (N+1) * (N+2) / 2;
	Real impPAv, admPAv, imp1Av, adm1Av;
	Real *rhsEx, *rhsEy, *rhsEz, *rhsHx, *rhsHy, *rhsHz;
	const Real *dEx, *dEy, *dEz, *dHx, *dHy, *dHz;
	const Real *dresEx, *dresEy, *dresEz, *dresHx, *dresHy, *dresHz;
	Real nx[np*nfp], ny[np*nfp], nz[np*nfp];
	Real rnx[np*nfp], rny[np*nfp], rnz[np*nfp];
	Real cnx[np*nfp], cny[np*nfp], cnz[np*nfp];
};

#endif /* SOLVERCURVEDFACE_H_ */
