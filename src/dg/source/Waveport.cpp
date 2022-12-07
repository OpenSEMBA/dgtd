#pragma once 

#include "Waveport.h"

namespace SEMBA::dgtd::dg::source {

////Math::Real
////DGWaveport::getNumericalGammaMGauss(
//// const Math::Real t,
//// const Math::Real dt,
//// const Math::Real amp,
//// const Math::Real delay,
//// const Math::Real spread,
//// const Math::Real kcm) const {
////	// Computes current step.
////	const size_t n = t / dt;
////	// Performs convolution.
////	size_t j;
////	Math::Real res = 0.0;
////	for (j = 0; j < n; j++) {
////		res +=
////		 getHm(j*dt, kcm)
////		 * getGauss((n-j)*dt, amp,delay,spread);
////	}
////	res *= dt;
////	res += getGaussDerivative(t, amp,delay,spread) / c0;
////	return res;
////}
//
//Math::Real DGWaveport::getHm(
// const Math::Real t,
// const Math::Real kcm) const {
//	if (t == 0) {
//		return (kcm * kcm * Constants::c0 / 2.0);
//	} else {
//		return (kcm / t * j1(kcm * Constants::c0 * t));
//	}
//}
//
//bool DGWaveport::checkNormalsAreEqual(
// const vector<pair<size_t, size_t> >& elemFace) const {
////	CVecR3 n1, n2;
////	for (size_t i = 1; i < elemFace.size(); i++) {
////		const size_t id1 = cells.getIdOfRelPos(elemFace[i-1].first);
////		const size_t f1 = elemFace[i-1].second;
////		n1 = cells.getPtrToCellWithId(id1)->getSideNormal(f1);
////		const size_t id2 = cells.getIdOfRelPos(elemFace[i].first);
////		const size_t f2 = elemFace[i].second;
////		n2 = cells.getPtrToCellWithId(id2)->getSideNormal(f2);
////		if (n1 != n2) {
////			return false;
////		}
////	}
////	return true;
//}

}