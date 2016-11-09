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
 * SolverPMLBiaxial.cpp
 *
 *  Created on: Aug 2, 2013
 *      Author: luis
 */

#include "../../dg/dispersive/DGPMLMultiaxial.h"

DGPMLMultiaxial::DGPMLMultiaxial(
        const PMVolumePML& mat,
        const CellGroup& cells,
        const bool useConductivity,
        const Real conductivity) : DGPML(mat, cells) {
//    J.set(dof, 0.0);
//    resJ.set(dof, 0.0);
//    rhsJ.set(dof, 0.0);
//    M.set(dof, 0.0);
//    resM.set(dof, 0.0);
//    rhsM.set(dof, 0.0);
}

DGPMLMultiaxial::~DGPMLMultiaxial() {

}

void DGPMLMultiaxial::addRHSToRes(
        const size_t e1, const size_t e2,
        const Real rka, const Real dt) {
//    size_t i, e;
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            resJ[i] *= rka;
//            resJ[i] += rhsJ[i] * dt;
//            resM[i] *= rka;
//            resM[i] += rhsM[i] * dt;
//        }
//    }
}

void DGPMLMultiaxial::updateWithRes(
        const size_t e1,
        const size_t e2,
        const Real rkb) {
//    size_t i, e;
//#ifdef SOLVER_USE_OPENMP
//#pragma omp parallel for private(i, e)
//#endif
//    for (i = 0; i < dof; i++) {
//        e = elem[(i / np) % nElem];
//        if (e1 <= e && e < e2) {
//            J[i] += resJ[i] * rkb;
//            M[i] += resM[i] * rkb;
//        }
//    }
}

