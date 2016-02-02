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
 * SolverDispersive.h
 *
 *  Created on: Sep 11, 2012
 *      Author: luis
 */

#ifndef SOLVERDISPERSIVEVOLUMIC_H_
#define SOLVERDISPERSIVEVOLUMIC_H_

using namespace std;

#include "DGDispersive.h"
#include "physicalModel/PMVolumeDispersive.h"

class DGDispersiveVolumic : public DGDispersive, public PMVolumeDispersive {
public:
    DGDispersiveVolumic();
    DGDispersiveVolumic(
            const PMVolumeDispersive&,
            const CellGroup& cells);
    virtual ~DGDispersiveVolumic();
    void computeRHSElectricPolarizationCurrents(
            const FieldR3& E,
            const UInt e1, const UInt e2);
    void computeRHSMagneticPolarizationCurrents(
            const FieldR3& H,
            const UInt e1, const UInt e2);
    void computeRHSElectric(
            FieldR3& rhsE,
            const FieldR3& E,
            const UInt e1, const UInt e2) const;
    void computeRHSMagnetic(
            FieldR3& rhsE,
            const FieldR3& E,
            const UInt e1, const UInt e2) const;
    void addRHSToRes(
            const UInt e1, const UInt e2,
            const Real rkb, const Real dt);
    void updateWithRes(
            const UInt e1,
            const UInt e2,
            const Real rkb);
    void addJumps(
            FieldR3& dE, FieldR3& dH,
            FieldR3& E, FieldR3& H,
            const UInt e1, const UInt e2);
private:
    static const UInt N = ORDER_N;
    static const UInt np = (N+1) * (N+2) * (N+3) / 6;
    //	PMVolumeDispersive mat;
    UInt nElem;
    UInt *elem;
    UInt dof, drudeDof;
    // Polarization currents. Size nK x Np x nPoles.
    // Data is stored in nK vectors of Np components for each pole.
    // First nK x Np data correspond to the first pole, and so on.
    FieldC3 P, J;
    FieldC3 rhsP, rhsJ;
    FieldC3 resP, resJ;
private:
    void build();
};
#endif /* SOLVERDISPERSIVE_H_ */
