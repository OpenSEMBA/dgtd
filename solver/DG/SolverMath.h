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
 * SolverMath.h
 *
 *  Created on: Sep 11, 2012
 *      Author: luis
 */
#ifndef SOLVERMATH_H_
#define SOLVERMATH_H_

// BLAS-L3 function for product. vecC = alpha * matA * vecB.
template<class T, int NR, int NC>
void am_v_prod(T* vC, const T *matA, const T *vB, const T alpha);
// BLAS-L3 function for product. vecC = -matA * vecB.
template<class T, int NR, int NC>
void minus_m_v_prod(T* vC, const T *matA, const T *vB);
// BLAS-L3 function for product. vecC = matA * vecB.
template<class T, int NR, int NC>
void m_v_prod(T* vC, const T *matA, const T *vB);
// BLAS-L3 function for product. vecC += alpha * matA * vecB.
template<class T, int NR, int NC>
void add_am_v_prod(T* vC, const T *matA, const T *vB, const T alpha);
// BLAS-L3 function for product. vecC += matA * vecB.
template<class T, int NR, int NC>
void add_m_v_prod(T* vC, const T *matA, const T *vB);
// BLAS-L3 function for product. vecC -= matA * vecB.
template<class T, int NR, int NC>
void sub_am_v_prod(T* vC, const T *matA, const T *vB, const T alpha);
// BLAS-L3 function for product. vecC -= matA * vecB.
template<class T, int NR, int NC>
void sub_m_v_prod(T* vC, const T *matA, const T *vB);
// BLAS-L3 function for product. vecC -= alpha * vecB.
template<class T, int NR>
void sub_a_v_prod(T* vC, const T *vB, const T alpha);

#include "SolverMath.hpp"

#endif

