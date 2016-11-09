/*
 * CellTri3.h
 *
 *  Created on: Jul 26, 2013
 *      Author: luis
 */

#ifndef CELLTRI3_H_
#define CELLTRI3_H_

#include "../../../common/geometry/Surface.h"
#include "../../dgtd/core/CellTri.h"

template<int TRI_N>
class CellTri3 : public CellTri<TRI_N>, public Tri3 {
public:
	static const size_t np = (TRI_N+1) * (TRI_N+2) / 2;
	static const size_t nfp = (TRI_N + 1);
	static const size_t faces = 3;
	static const size_t vertices = 3;
	CellTri3(
	 const Tri3& base);
	virtual
	~CellTri3();
	CVecC3
	 getRadiatedField(
	  const CVecC3 electric[np],
	  const CVecC3 magnetic[np],
	  const double frequency,
	  const pair<double, double> direction) const;
};

#include "../../dgtd/core/CellTri3.hpp"

#endif /* CELLTRI3_H_ */
