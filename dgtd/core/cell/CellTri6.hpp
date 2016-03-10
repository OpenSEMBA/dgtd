/*
 * CellTri6.cpp
 *
 *  Created on: Jul 26, 2013
 *      Author: luis
 */

#include "../../dgtd/core/CellTri6.h"

template<int TRI_N>
CellTri6<TRI_N>::CellTri6(
 const Tri6& base) {
	assert(this->tri.ncp == geo.ncp);
	id = base.getId();
	normal = base.getNormal();
	for (UInt i = 0; i < geo.np; i++) {
		v[i] = base.getV(i);
	}
}

template<int TRI_N>
CellTri6<TRI_N>::~CellTri6() {
	// TODO Auto-generated destructor stub
}

template<int TRI_N>
inline CVecC3
CellTri6<TRI_N>::getRadiatedField(
 const CartesianVector<complex<double>, 3> electric[np],
 const CartesianVector<complex<double>, 3> magnetic[np],
 const double frequency,
 const pair<double, double> direction) const {
	// Computes phase contribution.
	complex<double> phase[geo.ncp];
	const double beta = 2.0 * Constants::pi * frequency / SPEED_OF_LIGHT;
	complex<double> phaseShift(0.0, beta);
	SphericalVector sphDir(direction.first, direction.second);
	CVecR3 dir = sphDir.convertToCartesian();
	CVecR3 cNode[geo.ncp];
	getCubatureNodes(cNode);
	for (UInt c = 0; c < geo.ncp; c++) {
		phase[c] = exp(phaseShift * (double) dir.dot(cNode[c]));
	}
	// Computes integral.
	const double c0mu0 = SPEED_OF_LIGHT * VACUUM_PERMEABILITY;
	CVecR3 cNormal[geo.ncp];
	getCubatureNormals(cNormal);
	double csdf[geo.ncp];
	getCubatureDifferentials(csdf);
	CVecC3 res;
	CVecC3 complexDir;
	complexDir = MathUtils::convertToComplex(dir);
	CartesianVector<complex<double>, 3>	 complexNormal[geo.ncp];
	MathUtils::convertToComplex(complexNormal,cNormal,geo.ncp);
	static const SimplexTri<TRI_N> tri;
	CVecC3 M, J, integrand;
	for (UInt j = 0; j < tri.np; j++) {
		for (UInt c = 0; c < geo.ncp; c++) {
			M = - (complexNormal[c] ^ electric[j]);
			J =   (complexNormal[c] ^ magnetic[j]);
			integrand = ((J ^ complexDir) * c0mu0 + M) ^ complexDir;
			res += integrand *
			 (phase[c] * geo.cw[c] * csdf[c] / 2.0) * tri.ca[j][c];
		}
	}
	//res *= (phaseShift / (Constants::pi * 4.0)) * exp(phaseShift);
	return res;
}
