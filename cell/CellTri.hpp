/*
 * CellTri.cpp
 *
 *  Created on: Jul 26, 2013
 *      Author: luis
 */

#include "../../dgtd/core/CellTri.h"

template<int TRI_N>
CellTri<TRI_N>::CellTri() {
	// TODO Auto-generated constructor stub

}

template<int TRI_N>
CellTri<TRI_N>::~CellTri() {
	// TODO Auto-generated destructor stub
}

template<int TRI_N>
size_t
CellTri<TRI_N>::getNumberOfVertices() const {
	return vertices;
}

template<int TRI_N>
CVecR3
CellTri<TRI_N>::getSideNormal(const size_t s) const {
	cout << "ERROR @ CellTri:: Not Done" << endl;
	exit(-1);
}

