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

#include "SolverMath.h"

template<class T, int NR, int NC>
inline void am_v_prod(T* vC, const T *matA, const T *vB, const T alpha) {
	// Performs vC = alpha * matA * vB. Assumes row major ordering in mat.
	register int i, j, k = 0;
	for (i = 0; i < NR; i++) {
		vC[i] = (T) 0.0;
		for (j = 0; j < NC; j++) {
			vC[i] += matA[k++] * vB[j];
		}
		vC[i] *= alpha;
	}
}
 
template<class T, int NR, int NC>
inline void minus_m_v_prod(T* vC, const T *matA, const T *vB) {
	// Performs vecC = -matA * vecB. Assumes row major ordering in mat.
	register int i, j, k = 0;
	for (i = 0; i < NR; i++) {
		vC[i] = (T) 0.0;
		for (j = 0; j < NC; j++) {
			vC[i] -= matA[k++] * vB[j];
		}
	}
}
 
template<class T, int NR, int NC>
inline void m_v_prod(T* vC, const T *matA, const T *vB) {
	// Performs vC = alpha * matA * vB. Assumes row major ordering in mat.
	register int i, j, k = 0;
	for (i = 0; i < NR; i++) {
		vC[i] = 0.0;
		for (j = 0; j < NC; j++) {
			vC[i] += matA[k++] * vB[j];
		}
	}
}
 
template<class T, int NR, int NC>
inline void add_am_v_prod(T* vC, const T *matA, const T *vB, const T alpha) {
	// Performs vC += alpha * matA * vB. Assumes row major ordering in mat.
	register int i, j, k = 0;
	T sum;
	for (i = 0; i < NR; i++) {
		sum = (T) 0.0;
		for (j = 0; j < NC; j++) {
			sum += matA[k++] * vB[j];
		}
		vC[i] += alpha * sum;
	}
}
 
template<class T, int NR, int NC>
inline void add_m_v_prod(T* vC, const T *matA, const T *vB) {
	// Performs vC += matA * vB. Assumes row major ordering in mat.
	register int i, j, k = 0;
	for (i = 0; i < NR; i++) {
		for (j = 0; j < NC; j++) {
			vC[i] += matA[k++] * vB[j];
		}
	}
}

template<class T, int NR, int NC>
inline void sub_am_v_prod(T* vC, const T *matA, const T *vB, const T alpha) {
	// Performs vC += alpha * matA * vB.
	// Assumes row major ordering in mat.
	register int i, j, k = 0;
	T sum;
	for (i = 0; i < NR; i++) {
		sum = (T) 0.0;
		for (j = 0; j < NC; j++) {
			sum += matA[k++] * vB[j];
		}
		vC[i] -= alpha * sum;
	}
}
 
template<class T, int NR, int NC>
inline void sub_m_v_prod(T* vC, const T *matA, const T *vB) {
	// Performs vC -= alpha * matA * vB.
	// Assumes row major ordering in mat.
	register int i, j, k = 0;
	for (i = 0; i < NR; i++) {
		for (j = 0; j < NC; j++) {
			vC[i] -= matA[k++] * vB[j];
		}
	}
}

template<class T, int NR>
inline void sub_a_v_prod(T* vC, const T *vB, const T alpha) {
	for (int i = 0; i < NR; i++) {
		vC[i] -= alpha * vB[i];
	}
}
