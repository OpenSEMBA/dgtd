/*
 * CellTri6.h
 *
 *  Created on: Jul 26, 2013
 *      Author: luis
 */

#ifndef CELLTRI6_H_
#define CELLTRI6_H_

#include "../../../common/geometry/Surface.h"
#include "../../dgtd/core/CellTri.h"

template <int TRI_N>
class CellTri6 : public CellTri<TRI_N>, public Tri6 {
	static const UInt np = (TRI_N+1) * (TRI_N+2) / 2;
	static const UInt nfp = (TRI_N + 1);
	static const UInt faces = 3;
	static const UInt vertices = 3;
public:
	CellTri6(
	 const Tri6& base);
	virtual
	~CellTri6();
	CVecC3
	 getRadiatedField(
	  const CVecC3 electric[np],
	  const CVecC3 magnetic[np],
	  const double frequency,
	  const pair<double, double> direction) const;
};

#include "../../dgtd/core/CellTri6.hpp"

#endif /* CELLTRI6_H_ */
