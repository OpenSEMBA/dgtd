/*
 * Tri3.cpp
 *
 *  Created on: Jul 26, 2013
 *      Author: luis
 */

#include "Triangle3.h"

namespace SEMBA {
namespace Cudg3d {
namespace Cell {

template<int N>
Triangle3<N>::Triangle3(const Geometry::Tri3& base)
: Geometry::Tri3(base) {

}

template<int N>
Triangle3<N>::~Triangle3() {
	// TODO Auto-generated destructor stub
}
//
//template<int TRI_N>
//inline Math::CVecC3 Triangle3<TRI_N>::getRadiatedField(
// const Math::CVecC3 electric[np],
// const Math::CVecC3  magnetic[np],
// const double frequency,
// const std::pair<double, double> direction) const {
//	// Computes phase contribution.
//	complex<double> phase[geo.ncp];
//	const double beta = 2.0 * Math::Constants::pi * frequency / SPEED_OF_LIGHT;
//	complex<double> phaseShift(0.0, beta);
//	SphericalVector sphDir(direction.first, direction.second);
//	CVecR3 dir = sphDir.convertToCartesian();
//	CVecR3 cNode[geo.ncp];
//	getCubatureNodes(cNode);
//	for (size_t c = 0; c < geo.ncp; c++) {
//		phase[c] = exp(phaseShift * (double) dir.dot(cNode[c]));
//	}
//	// Computes integral.
//	const double c0mu0 = SPEED_OF_LIGHT * VACUUM_PERMEABILITY;
//	CVecR3 cNormal[geo.ncp];
//	getCubatureNormals(cNormal);
//	double csdf[geo.ncp];
//	getCubatureDifferentials(csdf);
//	CVecC3 res;
//	CVecC3 complexDir;
//	complexDir = MathUtils::convertToComplex(dir);
//	CartesianVector<complex<double>, 3>	 complexNormal[geo.ncp];
//	MathUtils::convertToComplex(complexNormal,cNormal,geo.ncp);
//	static const SimplexTri<TRI_N> tri;
//	CVecC3 M, J, integrand;
//	for (size_t j = 0; j < tri.np; j++) {
//		for (size_t c = 0; c < geo.ncp; c++) {
//			M = - (complexNormal[c] ^ electric[j]);
//			J =   (complexNormal[c] ^ magnetic[j]);
//			integrand = ((J ^ complexDir) * c0mu0 + M) ^ complexDir;
//			res += integrand *
//			 (phase[c] * geo.cw[c] * csdf[c] / 2.0) * tri.ca[j][c];
//		}
//	}
//	//res *= (phaseShift / (Constants::pi * 4.0)) * exp(phaseShift);
//	return res;
//}

}  /* namespace Cell */
}  /* namespace Cudg3d */
}  /* namespace SEMBA */


