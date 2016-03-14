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
 * SolverWaveport.cpp
 *
 *  Created on: Aug 26, 2013
 *      Author: luis
 */

#include "../../dg/sources/DGWaveport.h"

DGWaveport::DGWaveport() {
	posSF = NULL;
	posTFNB = NULL;
	posTF = NULL;
	gauss = NULL;
	hm = NULL;
}

DGWaveport::~DGWaveport() {
	// TODO Auto-generated destructor stub
}
//
//Real
//DGWaveport::getNumericalGammaMGauss(
// const Real t,
// const Real dt,
// const Real amp,
// const Real delay,
// const Real spread,
// const Real kcm) const {
//	// Computes current step.
//	const size_t n = t / dt;
//	// Performs convolution.
//	size_t j;
//	Real res = 0.0;
//	for (j = 0; j < n; j++) {
//		res +=
//		 getHm(j*dt, kcm)
//		 * getGauss((n-j)*dt, amp,delay,spread);
//	}
//	res *= dt;
//	res += getGaussDerivative(t, amp,delay,spread) / c0;
//	return res;
//}

Real DGWaveport::getHm(
 const Real t,
 const Real kcm) const {
	if (t == 0) {
		return (kcm * kcm * Constants::c0 / 2.0);
	} else {
		return (kcm / t * j1(kcm * Constants::c0 * t));
	}
}

bool DGWaveport::checkNormalsAreEqual(
 const vector<pair<size_t, size_t> >& elemFace) const {
//	CVecR3 n1, n2;
//	for (size_t i = 1; i < elemFace.size(); i++) {
//		const size_t id1 = cells.getIdOfRelPos(elemFace[i-1].first);
//		const size_t f1 = elemFace[i-1].second;
//		n1 = cells.getPtrToCellWithId(id1)->getSideNormal(f1);
//		const size_t id2 = cells.getIdOfRelPos(elemFace[i].first);
//		const size_t f2 = elemFace[i].second;
//		n2 = cells.getPtrToCellWithId(id2)->getSideNormal(f2);
//		if (n1 != n2) {
//			return false;
//		}
//	}
//	return true;
}
