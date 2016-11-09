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
 * SolverPMLUniaxial.h
 *
 *  Created on: Aug 2, 2013
 *      Author: luis
 */

#ifndef SOLVERPMLUNIAXIAL_H_
#define SOLVERPMLUNIAXIAL_H_

#include "../../dg/dispersive/DGPML.h"

template<Int D>
class DGPMLUniaxial : public DGPML {
public:
    DGPMLUniaxial(
            const PMVolumePML& mat,
            const CellGroup& cells,
            const bool useConductivity,
            const Real conductivity);
    virtual ~DGPMLUniaxial();
    void addRHSToRes(
            const size_t e1, const size_t e2,
            const Real rka, const Real dt);
    void updateWithRes(
            const size_t e1,
            const size_t e2,
            const Real rkb);
    void computeRHSElectric(
            FieldR3& rhsE,
            const FieldR3& E,
            const size_t e1, const size_t e2) const;
    void computeRHSMagnetic(
            FieldR3& rhsH,
            const FieldR3& H,
            const size_t e1, const size_t e2) const;
    void computeRHSElectricPolarizationCurrents(
            const FieldR3& E,
            const size_t e1, const size_t e2);
    void computeRHSMagneticPolarizationCurrents(
            const FieldR3& H,
            const size_t e1, const size_t e2);
protected:
    FieldR1 J, M, resJ, resM, rhsJ, rhsM;
private:
    bool check() const;
    static const CartesianAxis dir1 = CartesianAxis(D);
    static const CartesianAxis dir2 = CartesianAxis((D + 1)%3);
    static const CartesianAxis dir3 = CartesianAxis((D + 2)%3);
};

#include "../../dg/dispersive/DGPMLUniaxial.hpp"

typedef DGPMLUniaxial<x> DGPMLx;
typedef DGPMLUniaxial<y> DGPMLy;
typedef DGPMLUniaxial<z> DGPMLz;

#endif /* SOLVERPMLUNIAXIAL_H_ */
